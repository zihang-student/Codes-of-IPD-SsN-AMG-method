function [zeta,itamg,resamg,info] = Hybrid_AMG_detail(prob_data,amg_options)
% Original system:
%                   He * zeta = z
% with He = bk1*I + 1/tk*H, H = T + H0, and
%         S = diag(s), T = diag(t), H0 = A*S*A'.
% Let zeta = Q0 * u, we have an equivalent system:
%                     Ae * u = f
% with Ae = bk1*Q0^2 + 1/tk*(K+A0) and f = Q0*z, where
%         Q0 = diag(q,-p), K = Q0*T*Q0, A0 = Q0*H0*Q0.

%% He * zeta = z
bk1 = prob_data.bk1;tk = prob_data.tk;
q = prob_data.q;p = prob_data.p;
H0 = prob_data.H0;z = prob_data.z;
T = prob_data.T;M = size(H0,1);

amg_options.disp = [1,0,1];amg_options.time = 0;
%% Ae * u = f
qp = [q;-p];
if any(qp==0)
    error('p or q contains 0 !!!!!');
else
    Q0 = spdiags(qp,0,M,M);
end
A0 = Q0*H0*Q0; Q = Q0*Q0; K = Q0*T*Q0;
Ae = bk1*Q + 1/tk*(K+A0);f = Q0*z;
%% Check the connected components of H0/A0
% This is from the Lean-AMG toolbox
[blocks,sizes,ps,rs] = components(A0); num_comp = length(sizes);
%% Case 1
dK = diag(K);n = length(q);
if length(sizes) == 1 % A0/H0 is connected
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(dK)    % SPD
        % Notice that if sum(dK(pk)) > 0 but it is small then it is
        % still nearly singular.
        amg_options.isnsp = 0;
    else % nearly singular
        amg_options.isnsp = 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    amg_options.fnode = n;amg_options.guess = bk1*tk*rand(size(f));
    [u,itamg,resamg] = Class_AMG_detail(Ae,f,amg_options);
    if resamg > amg_options.retol
        amg_options.isnsp = ~amg_options.isnsp;
        [u,itamg,resamg] = Class_AMG_detail(Ae,f,amg_options);
    end
    it_num = 1;
end
%% Otherwise Case 2
if length(sizes) > 1  % A0/H0 is disconnected
    N0 = 1e2;itamg = 0;resamg = 0;it_num = 0;
    % ===================large components==========================
    large_size = sizes > N0;large_size = find(large_size);
    if ~isempty(large_size)
        for k = large_size
            pk = ps(rs(k):rs(k+1)-1);
            Aek = Ae(pk,pk);fk = f(pk);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if sum(dK(pk))   % SPD
                % Notice that if sum(dK(pk)) > 0 but it is small then it is
                % still nearly singular.
                amg_options.isnsp = 0;
            else % nearly singular
                amg_options.isnsp = 1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            amg_options.fnode = sum(pk<=n);
            amg_options.guess = bk1*tk*rand(size(fk));
            [dk,itamgk,resamgk] = Class_AMG_detail(Aek,fk,amg_options);
            if resamgk > amg_options.retol
                amg_options.isnsp = ~amg_options.isnsp;
                [dk,itamgk,resamgk] = Class_AMG_detail(Aek,fk,amg_options);
            end
            
            u(pk,1) = dk;itamg = max(itamg,itamgk);
            resamg = max(resamg,resamgk);
        end
        it_num = k;
    end
    % ===================large components==========================
    
    % ===================small components==========================
    blocks2sizes = sizes(blocks);small_comp = find(blocks2sizes<=N0);
    bb = blocks(small_comp);[~,b] = sort(bb);pk = small_comp(b);
    if ~isempty(pk)
        A0 = A0(pk,pk); Q = Q(pk,pk); K = K(pk,pk);
        f = f(pk);Ae = bk1*Q + 1/tk*(K+A0);
        % Method I: Direct Is this robust ???
        u(pk,1) = Ae \ f;
        % disp(norm(u1(pk,1)));
        %===========================
        % Method II: PCG Does this work well ??
        %         Y = sparse(1:M,blocks,ones(M,1),M,num_comp);
        %         Y = Y(pk,bb);QK = bk1*Q + 1/tk*K; % A0*Y = 0
        %         augAe = [(Y')*QK*Y Y'*QK;QK*Y Ae];augf = [Y'*f;f];
        %
        %         pcg_options.guess = [zeros(size(Y,2),1);zeros(size(f))];
        %         pcg_options.retol = 1e-11;pcg_options.maxit = 1e4;
        %         pcg_options.precd = 2;
        %
        %         [U,itpcg,respcg] = PCG(augAe,augf,pcg_options);
        %         u(pk,1) = Y*U(1:size(Y,2)) + U(size(Y,2)+1:end);
        %         itamg = max(itamg,itpcg);resamg = max(resamg,respcg);
    end
    % ===================small components==========================
end
%% zeta
zeta = Q0*u; info = [num_comp,it_num];
end