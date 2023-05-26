function [x,it,rel_res,rel_resk,rhok] = twogrid(A,b,amg_options)
%% Two-grid method for solving the linear algebraic system
%                        Ax = b, A = eps*I + K + A0,
% where K>=0 is diagonal and A0 is a connected graph Laplacian.

% amg_options
% .retol: relative tolerance
% .bigph: A0 is a bigraph (1) or not (0)
% .maxit: maximal number of iterations
% .smoth: number of smoothing iterations
% .isnsp: is nearly singular problem
% .guess: initial guess x0

% .fnode: number of fine nodes (This is only for bigraph)

if nargin == 2
    % default setting------------------------------------------------------
    amg_options = struct('retol',1e-12,'bigph',0,'maxit',2e1,...
        'smoth',10,'isnsp',1,'guess',zeros(size(b)));
    % default setting------------------------------------------------------
end
if isempty(amg_options.retol);amg_options.retol = 0;end
if isempty(amg_options.bigph);amg_options.bigph = 0;end
if isempty(amg_options.maxit);amg_options.maxit = 5e1;end
if isempty(amg_options.smoth);amg_options.smoth = 3;end
if isempty(amg_options.isnsp);amg_options.isnsp = 0;end
if isempty(amg_options.guess);amg_options.guess = zeros(size(b));end

if isfield(amg_options,'fnode')
    if isempty(amg_options.fnode)
        amg_options.fnode = 0;end
else
    amg_options.fnode = 0;
end

if amg_options.bigph && (amg_options.fnode <= 0)
    error('bigph = 1 requires fnode > 0');
end

%% Setup Phase
N = size(A,1);
% 1. smoother R (inverse)
if amg_options.bigph  % A0 is a bigraph
    Nf = amg_options.fnode;Nc = N - Nf;
    Aff = A(1:Nf,1:Nf);Afc = A(1:Nf,Nf+1:N);
    if norm(Aff-spdiags(diag(Aff),0,Nf,Nf),'inf')
        error('Nf is not right for the bigraph');
    end
    Acc = A(Nf+1:N,Nf+1:N);
    invV = spdiags(1./diag(Aff),0,Nf,Nf);
    invT = spdiags(1./diag(Acc),0,Nc,Nc);
    % =======================GS=====================
    O = spalloc(Nf,Nc,1);R = [invV O;-invT*Afc'*invV invT];  % R^{-1} !!!!!
    % =======================GS=====================
    % ======================SSOR====================
    %     R = w*(2-w)*[invV + w^2*invV*Afc*invT*Afc'*invV,...
    %         -w*invV*Afc*invT;-w*invT*Afc'*invV, invT]; % R^{-1} !!!!!
    % ======================SSOR====================
else
    % ================weighted Jacobi===============
    R = .5*spdiags(1./diag(A),0,N,N); % R^{-1} !!!!!
    % ================weighted Jacobi===============
end
% 2. prolongation and coarse matrix
if amg_options.bigph % A0 is a bigraph
    % ideal interpolation
    W = - Aff\Afc;
    if amg_options.isnsp == 1
        D = spdiags(W*ones(Nc,1),0,Nf,Nf);W = D\W;
    end
    Pro = [W;speye(Nc)];
else
    [indC,indF,As] = mis_set(A,1/4);
    %% Transfer operators
    C_node = find(indC);Nc = size(C_node,1);
    F_node = find(indF);Nf = size(F_node,1);
    %%%%%%%===============NEW VERSION===========%%%%%%%%
    p = [F_node;C_node];AA = A(p,p);
    Aff = AA(1:Nf,1:Nf);Afc = AA(1:Nf,Nf+1:N);
    
    Dff = spdiags(diag(Aff),0,Nf,Nf);W1 = - Dff\Afc;
    as = speye(Nf) + As(F_node,F_node);
    Affs = Aff.*as;W2 = - Dff\Affs*W1;
    W = W1 + W2;
    % Be careful with the case W = 0 !!!!!!!
    if ~isempty(W*ones(Nc,1) == 0)
        W = W1 + 0.5*W2;
    end
    if amg_options.isnsp == 1
        D = spdiags(W*ones(Nc,1),0,Nf,Nf);W = D\W;
    end
    P = [W;speye(Nc)];Pro = P;Pro(p,:) = P;
    %%%%%%%===============NEW VERSION===========%%%%%%%%
end
Ac = Pro'*A*Pro;

%% Solve Phase
it = 0;rhok = nan(amg_options.maxit,1);
rel_resk = ones(amg_options.maxit,1);
x = amg_options.guess;res0 = norm(A*x-b);

if res0 == 0
    rel_res = 0;rel_resk = 0;rhok = inf;
else
    it = it + 1;
    while rel_resk(it) > amg_options.retol && it <= amg_options.maxit
        r = b - A*x;
        x = x + twogrid_iteration(A,r,Ac,R,Pro,amg_options);
        res = norm(A*x-b);rel_res = res/res0;
        rel_resk(it+1) = rel_res;
        rhok(it+1) = res/norm(r);
        it = it + 1; if rhok(it) > 1; break;end
    end
    rel_resk = rel_resk(1:it);rhok = rhok(1:it);it = it - 1;
end
end

function [e] = twogrid_iteration(A,r,Ac,R,Pro,amg_options)
N = size(A,1);Rt = R';
isnsp = amg_options.isnsp;
smoth = amg_options.smoth;
% Presmoothing
e = zeros(size(r));
if  isnsp % add kernel space
    xi = ones(N,1);Axi = A*xi;xx = xi'*Axi;
    for i = 1:smoth
        g = r-A*e;xig = xi'*g;
        g = xi*(xig/xx) + R*(g-Axi*(xig/xx));
        e = e + g;
    end
else
    for i = 1:smoth; e = e + R*(r-A*e); end
end
% Restriction
rr = r-A*e;rrc = Pro'*rr;
% Correction % Ac * eec = rrc
pcg_options = struct('retol',[],'maxit',1e2,'precd',2,'guess',[]);
[eec] = PCG(Ac,rrc,pcg_options);
% Prolongation
ee = Pro*eec;e = e + ee;
% Postmoothing
if isnsp  % add kernel space
    for i = 1:smoth
        g = r-A*e;xig = xi'*g;
        g = xi*(xig/xx) + Rt*(g-Axi*(xig/xx));
        e = e + g;
    end
else
    for i = 1:smoth;e = e + Rt*(r-A*e); end
end
end