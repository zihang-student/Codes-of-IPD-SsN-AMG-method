

% Original system:
%                   He * zeta = z
% with He = bk1*I + 1/tk*H, H = T + H0, and
%         S = diag(s), T = diag(t), H0 = A*S*A'.

add_path;

prob = 2;m = 4e3;n = m; % Optimal mass transport
[c,r,l,p,q,gama] = prob_set(m,n,prob);
b = [r;l];nb = norm(b);
A = [kron(speye(n),p');kron(q',speye(m))];

bk1 = 1e-7;tk = 1e-2;
t = zeros(m+n,1);
% sp = 0.005;s = randn(m*n,1);s(abs(s)<sp) = abs(s(abs(s)<sp));s(abs(s)>=sp) = 0;
Sinkhorn;s = X(:);

% S = reshape(s,n,m);
% S(:,[1,3,6,10,20,50]) = 0;
% S(:,100) = 0;S(200,100) = 1;
% % S(200,:) = 0;S(200,100) = 1;
% s = S(:);

T = spdiags(t,0,m+n,m+n);
S = spdiags(s,0,m*n,m*n);


H0 = A*S*A';H = T + H0;
He = bk1*speye(m+n) + 1/tk*H;
z = He*He(:,1);
%% Problem data
prob_data.bk1 = bk1;prob_data.tk = tk;
prob_data.q = q;prob_data.p = p;prob_data.z = z;
prob_data.H0 = H0;prob_data.T = T;
% save strange_data prob_data;
%% Hybrid_AMG/Twogrid
% 1. Twogrid
amg_options = struct('retol',1e-11,'bigph',1,'maxit',5e1,'theta',1/4,...
    'smoth',5,'cycle','v','isnsp',1,'inter',1,'guess',[]);

tic;[~,itamg,resamg] = Hybrid_twogrid(prob_data,amg_options);
disp(['Hybrid_twogrid: it = ', num2str(itamg,'%2d'),'  rel_res = ',...
    num2str(resamg,'%4.2e'),'  time used ',num2str(toc,'%02.2f'),'s']);
% 2. AMG
tic;[zeta,itamg,resamg] = Hybrid_AMG(prob_data,amg_options);
disp(['Hybrid_AMG:     it = ', num2str(itamg,'%2d'),'  rel_res = ',...
    num2str(resamg,'%4.2e'),'  time used ',num2str(toc,'%02.2f'),'s']);
% 3. AMG with more details
amg_options.disp = [1,0,1];amg_options.time = 0;
disp('Hybrid_AMG_detail====================================================================');
[~,itamg,resamg] = Hybrid_AMG_detail(prob_data,amg_options);
disp('Hybrid_AMG_detail====================================================================');
%% PCG
pcg_options = struct('retol',1e-11,'maxit',1e4,'precd',2,'guess',[]);
tic;[d1,it,res] = PCG(He,z,pcg_options);
disp(['Precond_CG: it = ', num2str(it,'%2d'),'  rel_res = ',num2str(res,'%4.2e'),...
    '  time used ',num2str(toc,'%02.2f'),'s']);
% %% Hybrid_PCG
% tic;[~,itpcg,respcg,info] = Hybrid_PCG(prob_data,pcg_options);
% disp(['Hybrid_PCG: it = ', num2str(itpcg,'%2d'),'  rel_res = ',...
%     num2str(respcg,'%4.2e'),'  time used ',num2str(toc,'%02.2f'),'s']);
%% Aug_PCG
tic;[d2,itpcg,respcg,info] = aug_PCG(prob_data,pcg_options);
disp(['Augmnt_PCG: it = ', num2str(itpcg,'%2d'),'  rel_res = ',...
    num2str(respcg,'%4.2e'),'  time used ',num2str(toc,'%02.2f'),'s']);
%%
disp('    _______________________________________________');
disp('     #component   #small(direct)    #large(iter)   ');
disp('    -----------------------------------------------');
aa = '       %4d          %4d              %4d\n';
fprintf(aa,info);
disp('    -----------------------------------------------');