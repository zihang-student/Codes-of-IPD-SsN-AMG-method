function [x,it,rel_res,rel_resk,rhok] = Class_AMG(A,b,amg_options)
%% Classical AMG for solving the linear algebraic system
%                        Ax = b, A = eps*I + K + A0,
% where K>=0 is diagonal and A0 is a connected graph Laplacian.

% amg_options
% .retol: relative tolerance
% .bigph: A0 is a bigraph or not
% .maxit: maximal number of iterations
% .theta: strength connection parameter
% .smoth: number of smoothing iterations
% .cycle: 'w' for W-cycle or 'v' for V-cycle
% .isnsp: is nearly singular problem
% .guess: initial guess x0
% .inter: interpolation choice
%         w = 0 direct interpolation
%         w = 1 standard interpolation  % default choice
%         w = 2 ideal interpolation

if nargin == 2
    % default setting------------------------------------------------------
    amg_options = struct('retol',1e-12,'bigph',0,'maxit',2e1,'theta',1/4,...
        'smoth',10,'cycle',1,'isnsp',1,'inter',1,'guess',zeros(size(b)));
    % default setting------------------------------------------------------
end
if isempty(amg_options.retol); amg_options.retol = 1e-12;end
if isempty(amg_options.bigph); amg_options.bigph = 0;end
if isempty(amg_options.maxit); amg_options.maxit = 5e1;end
if isempty(amg_options.theta); amg_options.theta = 1/4;end
if isempty(amg_options.smoth); amg_options.smoth = 3;end
if isempty(amg_options.cycle); amg_options.cycle = 'v';end
if isempty(amg_options.isnsp); amg_options.isnsp = 0;end
if isempty(amg_options.inter); amg_options.inter = 1;end
if isempty(amg_options.guess); amg_options.guess = zeros(size(b));end

if amg_options.bigph
    if isempty(amg_options.fnode) || amg_options.fnode <= 0
        error('amg_options.bigph = 1 requires Nf > 0');
    end
end
%% Setup phase
clear  Ack Prok J smoth_it Rk;
global Ack Prok J smoth_it Rk;
smoth_it = amg_options.smoth;

J = 1;Ak = A;Pro = 1;
Ack{J} = Ak;Prok{J} = Pro;dofk = size(A,1);
if amg_options.bigph  % A0 is a bigraph
    % R{1} uses GS
    Nf = amg_options.fnode;Nc = dofk - Nf;
    V = Ak(1:Nf,1:Nf);U = Ak(1:Nf,Nf+1:dofk);
%     if norm(V-spdiags(diag(V),0,Nf,Nf),'inf')
%         error('Nf is not true for bigraph');
%     end
    T = Ak(Nf+1:dofk,Nf+1:dofk);
    invV = spdiags(1./diag(V),0,Nf,Nf);
    invT = spdiags(1./diag(T),0,Nc,Nc);
    O = spalloc(Nf,Nc,1);
    Rk{J} = [invV O;-invT*U'*invV invT];  % R^{-1} !!!!!
    % R{1} uses SSOR
    %     w = 1.6;Nf = amg_options.fnode;Nc = dofk - Nf;
    %     V = Ak(1:Nf,1:Nf);U = Ak(1:Nf,Nf+1:dofk);
    %     if norm(V-spdiags(diag(V),0,Nf,Nf),'inf')
    %         error('Nf is not true for bigraph');
    %     end
    %     T = Ak(Nf+1:dofk,Nf+1:dofk);
    %     invV = spdiags(1./diag(V),0,Nf,Nf);
    %     invT = spdiags(1./diag(T),0,Nc,Nc);
    %     Rk{J} = w*(2-w)*[invV + w^2*invV*U*invT*U'*invV,...
    %         -w*invV*U*invT;-w*invT*U'*invV, invT]; % R^{-1} !!!!!
else % Otherwise use (weighted) Jacobi
    Rk{J} = 0.5*spdiags(1./diag(Ak),0,size(Ak,1),size(Ak,1)); % R^{-1} !!!!!
end

% while size(Ak,1) > 1 + fix(size(A,1)/2)
while size(Ak,1) > 1 + fix(size(A,1)^(1/3))
% while size(Ak,1) > 5e2
    [Ak,Pro] = transfer(Ack{J},amg_options);
%     if max(max(isnan(Ak)))
%         disp(J);    error('Ak contains nan');
%     end
    J = J + 1; Ack{J} = Ak; Prok{J} = Pro;
    % R{J} uses (weighted) Jacobi
    Rk{J} = 0.5*spdiags(1./diag(Ak),0,size(Ak,1),size(Ak,1)); % R^{-1} !!!!!
end
%% Solve phase
it = 0;rhok = nan(amg_options.maxit,1);
rel_resk = ones(amg_options.maxit,1);
x = amg_options.guess;res0 = norm(A*x-b);

if res0 == 0
    rel_res = 0;rel_resk = 0;rhok = inf;
else
    it = it + 1;
    while rel_resk(it) > amg_options.retol && it <= amg_options.maxit
        r = b - A*x;
        if strcmp(amg_options.cycle,'v') % vcycle
            x = x + MG_Vcycle(r,amg_options.isnsp);
        end
        if strcmp(amg_options.cycle,'w') % wcycle
            x = x + MG_Wcycle(r,amg_options.isnsp);
        end
        res = norm(A*x-b);rel_res = res/res0;
        rel_resk(it+1) = rel_res;
        rhok(it+1) = res/norm(r);
        it = it + 1; if rhok(it) > 1; break;end
    end
    rel_resk = rel_resk(1:it);rhok = rhok(1:it);it = it - 1;
end
clear Ack Prok J smoth_it Rk;
end