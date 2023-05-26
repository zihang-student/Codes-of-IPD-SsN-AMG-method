function [x,it,rel_res,rel_resk,rhok] = twogrid_bigph(A,b,amg_options)
%% Two-grid method for solving the linear algebraic system
%                        Ax = b, A = eps*I + K + A0,
% where K>=0 is diagonal and A0 is a connected bigraph Laplacian.

% amg_options
% .retol: relative tolerance
% .maxit: maximal number of iterations
% .smoth: number of smoothing iterations
% .isnsp: is nearly singular problem
% .guess: initial guess x0

if nargin == 2
    % default setting------------------------------------------------------
    amg_options = struct('retol',1e-12,'maxit',2e1,'fnode',0,...
        'smoth',10,'isnsp',1,'guess',zeros(size(b)));
    % default setting------------------------------------------------------
end
if isempty(amg_options.retol);amg_options.retol = 0;end
if isempty(amg_options.maxit);amg_options.maxit = 5e1;end
if isempty(amg_options.smoth);amg_options.smoth = 3;end
if isempty(amg_options.isnsp);amg_options.isnsp = 0;end
if isempty(amg_options.guess);amg_options.guess = zeros(size(b));end

%% Setup Phase
N = size(A,1);
% 1. smoother R (inverse)
Nf = amg_options.fnode;Nc = N - Nf;
Aff = A(1:Nf,1:Nf);Afc = A(1:Nf,Nf+1:N);
Acc = A(Nf+1:N,Nf+1:N);
invV = spdiags(1./diag(Aff),0,Nf,Nf);
invT = spdiags(1./diag(Acc),0,Nc,Nc);
% =======================GS=====================
O = spalloc(Nf,Nc,1);R = [invV O;-invT*Afc'*invV invT];  % R^{-1} !!!!!
% =======================GS=====================
% ======================SSOR====================
%     w = 1.5;R = w*(2-w)*[invV + w^2*invV*Afc*invT*Afc'*invV,...
%         -w*invV*Afc*invT;-w*invT*Afc'*invV, invT]; % R^{-1} !!!!!
% ======================SSOR====================

% 2. prolongation and coarse matrix
% ideal interpolation
W = - Aff\Afc;
if amg_options.isnsp == 1
    D = spdiags(W*ones(Nc,1),0,Nf,Nf);W = D\W;
end
Pro = [W;speye(Nc)];

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
        x = x + twogrid_it(A,r,Ac,R,Pro,amg_options);
        res = norm(A*x-b);rel_res = res/res0;
        rel_resk(it+1) = rel_res;
        rhok(it+1) = res/norm(r);
        it = it + 1; if rhok(it) > 1; break;end
    end
    rel_resk = rel_resk(1:it);rhok = rhok(1:it);it = it - 1;
end
end

function [e] = twogrid_it(A,r,Ac,R,Pro,amg_options)
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
% 1. augmented PCG===================
% xxi = ones(size(Ac,1),1);
% augAc = [(xxi')*Ac*xxi xxi'*Ac;Ac*xxi Ac];augf = [xxi'*rrc;rrc];
% pcg_options = struct('retol',[],'maxit',5e1,'precd',2,'guess',[]);
% pcg_options.guess = [0;zeros(size(rrc))];
% [U] = PCG(augAc,augf,pcg_options);eec = xxi*U(1) + U(2:end);
% augmented PCG======================
% 2. PCG=============================
pcg_options = struct('retol',[],'maxit',1e2,'precd',2,'guess',[]);
[eec] = PCG(Ac,rrc,pcg_options);
% PCG=============================
% 3. pcg==========================
% M = spdiags(diag(Ac),0,size(Ac,1),size(Ac,1));
% eec = pcg(Ac,rrc,1e-10,1e2,M);
% pcg=============================
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