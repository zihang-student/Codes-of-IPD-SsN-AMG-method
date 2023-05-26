function [x,it,rel_res,rel_resk,rhok,ratio,level] = Class_AMG_detail(A,b,amg_options)
%% Classical AMG for solving the linear algebraic system
%                        Ax = b, A = eps*I + K + A0,
% where K>=0 is diagonal and A0 is a connected graph Laplacian. This
% version provides more details.

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

% .time:  record or not
% .disp:  display setup phase, coarse mesh or iteration

if nargin == 2
    % default setting------------------------------------------------------
    amg_options = struct('retol',1e-12,'bigph',0,'maxit',2e1,'theta',1/4,'smoth',10,...
        'cycle',1,'isnsp',1,'inter',1,'guess',zeros(size(b)),'disp',[0,0,0],'time',0);
    % default setting------------------------------------------------------
end
if isempty(amg_options.retol); amg_options.retol = 1e-12;end
if isempty(amg_options.bigph); amg_options.bigph = 0;end
if isempty(amg_options.maxit); amg_options.maxit = 5e1;end
if isempty(amg_options.theta); amg_options.theta = 1/40;end
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

if isempty(amg_options.time) ; amg_options.time  = 0;end
if isempty(amg_options.disp) ; amg_options.disp  = [0,0,0];end

if amg_options.time; amg_options.disp = [0,0,0]; end

flag_setup = amg_options.disp(1);
flag_mesh = amg_options.disp(2);
flag_iter = amg_options.disp(3);
%% Setup phase
clear  Ack Prok J smoth_it Rk;
global Ack Prok J smoth_it Rk;
smoth_it = amg_options.smoth;

J = 1;Ak = A;Pro = 1;
Ack{J} = Ak;Prok{J} = Pro;dofk = size(A,1);
if amg_options.bigph  % A0 is a bigraph
    % R{1} uses GS
    %     Rk{J} = tril(Ak); % R^{-1} !!!!!
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
    %         -w*invV*U*invT;-w*invT*U'*invV, invT];  % R^{-1} !!!!!
else % Otherwise use (weighted) Jacobi
    Rk{J} = 0.5*spdiags(1./diag(Ak),0,size(Ak,1),size(Ak,1)); % R^{-1} !!!!!
end

if flag_setup  % prepare setup phase details
    nnzAk = nnz(Ak);nnz_total = nnzAk;
    Ak0 = Ak - spdiags(diag(Ak),0,dofk,dofk);
    edge = sign(triu(Ak0));num_edge = nnz(edge);
    neg_edge = nnz(edge+abs(edge));
    deg_node = full(sum(Ak0~=0,2));
    
    disp('    __________________________________Setup_Phase____________________________________');
    disp('     level   #nnz    spa.    #node     #edge   #neg_edge   min_deg  ave_deg  max_deg   ');
    disp('    ---------------------------------------------------------------------------------');
    info1 = '       %2d    %6d   %02.2f    %5d     %6d      %2d        %3d       %2.0f      %3d\n';
    
    bb = [1,nnzAk,nnzAk/dofk^2,dofk,num_edge,neg_edge];
    aa = [bb,min(deg_node),sum(deg_node)/dofk,max(deg_node)];
    fprintf(info1,aa);tic;
end

if flag_mesh
    Sa{J} = [];indC{J} = [];
end

% while size(Ak,1) > 1 + fix(size(A,1)/2)
while size(Ak,1) > 1 + fix(size(A,1)^(1/3))
    % while size(Ak,1) > 2e2
    [Ak,Pro,Ask,indCk] = transfer(Ack{J},amg_options);
    J = J + 1; Ack{J} = Ak; Prok{J} = Pro;dofk = size(Ak,1);
    % R{J} uses (weighted) Jacobi
    Rk{J} = 0.5*spdiags(1./diag(Ak),0,size(Ak,1),size(Ak,1)); % R !!!!!
    
    if flag_setup % prepare setup phase details
        nnzAk = nnz(Ak);nnz_total  = nnz_total + nnzAk;
        Ak0 = Ak - spdiags(diag(Ak),0,dofk,dofk);
        edge = sign(triu(Ak0));num_edge = nnz(edge);
        neg_edge = nnz(edge+abs(edge));
        deg_node = full(sum(Ak0~=0,2));
        
        bb = [J,nnzAk,nnzAk/dofk^2,dofk,num_edge,neg_edge];
        aa = [bb,min(deg_node),sum(deg_node)/dofk,max(deg_node)];
        fprintf(info1,aa);
    end
    if flag_mesh  % prepare mesh details
        Sa{J-1} = Ask;indC{J} = indCk;
    end
end
if flag_setup
    disp('    ---------------------------------------------------------------------------------');
    ratio = nnz_total/nnz(A);
    fprintf('     total #nnz = %6d, ratio = %02.2f, time used %02.2fs\n',[nnz_total,ratio,toc]);
    disp('    ---------------------------------------------------------------------------------');
end
%% Show AMG mesh
if flag_mesh
    Ga = Ack;
    %     show_amg_mesh(Ga,Sa,indC,true); % fine(with strength connection) to coarse
    show_amg_mesh(Ga,Sa,indC); % fine to coarse
end
%% Solve phase
it = 0;rhok = nan(amg_options.maxit,1);
rel_resk = ones(amg_options.maxit,1);
x = amg_options.guess;res0 = norm(A*x-b);

if res0 == 0
    rel_res = 0;rel_resk = 0;rhok = inf;
else
    if flag_iter
        fprintf('\n');
        if strcmp(amg_options.cycle,'v') % vcycle
            disp('    _______AMG-Vcycle______');
        else
            disp('    _______AMG-Wcycle______');
        end
        disp('     iter   rel_res   c.f.');
        disp('    -----------------------');
        info3 = '     %2d    %4.2e   %02.2f\n';
        aa = [it,rel_resk(it+1),rhok(it+1)];
        fprintf(info3,aa);tic;
    end
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
        rel_resk(it+1) = rel_res;rhok(it+1) = res/norm(r);
        if flag_iter
            if res == 0
                aa = [it,rel_res,inf];fprintf(info3,aa);
            else
                aa = [it,rel_res,rhok(it+1)];fprintf(info3,aa);
            end
        end
        it = it + 1; if rhok(it) > 1; break;end
    end
    rel_resk = rel_resk(1:it);rhok = rhok(1:it);it = it - 1;
    if flag_iter
        disp('    -----------------------');
        fprintf('      AMG time used %02.2fs\n',toc);
        disp('    -----------------------');
    end
end
level = J; clear  Ack Prok J smoth_it Rk;
end