function [Ac,Pro,As,indC] = transfer(A,amg_options)
% Given a SYMMETRIC matrix A, this returns the transfer data from
% A to its coarse level Ac with classical interpolation methods:
%       w = 0 direct interpolation
%       w = 1 standard interpolation  % default choice
%       w = 2  ideal   interpolation

if nargin == 1
    amg_options = struct('theta',1/40,'bigph',0,'inter',1,'diag',0);
end

if isempty(amg_options.theta); amg_options.theta = 1/4;end
if isempty(amg_options.bigph); amg_options.bigph = 0;end
if isempty(amg_options.inter); amg_options.inter = 1;end
if isempty(amg_options.isnsp);  amg_options.isnsp = 0;end

global J; N = size(A,1);

if J == 1 && amg_options.bigph % A0 is a bigraph
    Nf = amg_options.fnode;Nc = N - Nf;
    Aff = A(1:Nf,1:Nf);Afc = A(1:Nf,Nf+1:N);W = -Aff\Afc;
    if amg_options.isnsp == 1
        D = spdiags(W*ones(Nc,1),0,Nf,Nf);W = D\W;
    end
    Pro = [W;speye(Nc)];
    if nargout > 2
        indC = false(N,1);indC(Nf+1:N) = true;
        As = strength(A) >= amg_options.theta;
    end
else
    % C/F splitting
    %         % Method I
    %         As = strength(A) >= amg_options.theta;
    %         [indC,indF] = cf_split(As); % logical
    %         if sum(indF) == 0
    %             N0 = min(fix(length(indF)/2)+1,10);
    %             dd = ceil(rand(N0,1)*length(indF)); % randomly chose N0 nodes
    %             indC = false(size(indC));indC(dd) = true; indF = ~indC;
    %         end
    % Method II
    [indC,indF,As] = mis_set(A,amg_options.theta);
    %% Transfer operators
    C_node = find(indC);Nc = size(C_node,1);
    F_node = find(indF);Nf = size(F_node,1);
    %%%%%%%===============NEW VERSION===========%%%%%%%%
    p = [F_node;C_node];AA = A(p,p);
    Aff = AA(1:Nf,1:Nf);Afc = AA(1:Nf,Nf+1:N);
    if amg_options.inter < 2
        Dff = spdiags(diag(Aff),0,Nf,Nf);W1 = - Dff\Afc;
        as = speye(Nf) + As(F_node,F_node);
        Affs = Aff.*as;W2 = - Dff\Affs*W1;
        W = W1 + amg_options.inter*W2;
        % Be careful with the case W = 0 !!!!!!!
        if ~isempty(W*ones(Nc,1) == 0)
            W = W1 + 0.5*W2;
        end
    else % ideal interpolation
        W = - Aff \ Afc;
    end
    if amg_options.isnsp == 1
        D = spdiags(W*ones(Nc,1),0,Nf,Nf);W = D\W;
    end
    P = [W;speye(Nc)];Pro = P;Pro(p,:) = P;
    %%%%%%%===============NEW VERSION===========%%%%%%%%
end
Ac = Pro'*A*Pro;
end