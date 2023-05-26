%% Accelerated ADMM for problem CLASS 1: Transport-like problem
%        min_{0<= x <=gama}  <c,x>  s.t. A*x = b,
% where gama = vec(Gama), c in R^{mn}_+, r \in R^n_+, l in R^m_+ and
%      [ In  ox 1m']        [ r ]
%  A = [           ]    b = [   ]
%      [ 1n' ox Im ]        [ l ]
% with the conditions: <r,1n> >= <l,1m>,l <= Gama*1n, r <= Gama'*1m.

%  1. Assignment problem: gama = inf, r = 1n, l = 1m
%  2. Optimal transport : gama = inf,
%  3. Capacity constrained mass transport: the general case

addpath '..';
%% 0. Problem Setting
% m = 1e2;n = m; prob = 0;       % half l1-norm = 0.5*norm(r-l,1)
% disp(['  0.5*||l-r||_1 = ', num2str(0.5*norm(l-r,1),'%4.5e')]);
% m = 5e2;n = m; prob = 1;  % Assignment problem
m = 50e2;n = m; prob = 2;  % Optimal transport
% m = 10e2;n = m; prob = 3;  % Capacity constrained mass transport

% [c,r,l,p,q,gama] = prob_set(m,n,prob);

% load('./InputData/data1-5000.mat');
% save('./InputData/data1-cap-1000.mat','p','q','r','l','m','n','c','gama')
% save('./InputData/data1-5000.mat','p','q','r','l','m','n','c','gama')


load('./InputData/apd_class1_l2_apd.mat');
APDtype=4;

m=apd_data(APDtype,1); n=m;
q = ones(n,1); p = ones(m,1);
Gama = + inf(m,n);gama = Gama(:);
c=cost_data{APDtype};
l=l_data{APDtype};
r=r_data{APDtype};



%% 1. Algorithm Setting
% 1.1 auxiliary vars
prox = @(x) min(max(0,x),gama);
b = [r;l];Atb = Aty(b,p,q);z0 = zeros(m+n,1);
% 1.2 parameters
muf = 0; gk = 1;bk = 1; maxit = 5e3; KKT_Tol = 1e-6;
% 1.3 initial guess
%%%%%%%%%% MATLAB linprog %%%%%%%%%%
% % options = optimoptions('linprog','Algorithm','interior-point-legacy',...
% %     'Display','iter','MaxIterations',100);
% options = optimoptions('linprog','Algorithm','interior-point-legacy',...
%     'Display','iter','MaxIterations',10,'OptimalityTolerance',1e-1);
% A = [kron(speye(n),p');kron(q',speye(m))];
% [mat_xk,fk,~,~,mat_lk] = linprog(c,[],[],A,b,spalloc(m*n,1,0),gama,[],options);
% mat_lk = mat_lk.eqlin;xk = mat_xk; vk = xk;wk = xk;pik = wk;lk = [mat_lk;xk];
% disp(['  fxk = ', num2str(fk,'%4.5e')]);
%%%%%%%%%% random initial %%%%%%%%%%
xk = spalloc(m*n,1,0);vk = xk;wk = xk;pik = wk;lk = [z0;xk];
% 1.4 residuals
fxk = - inf*ones(maxit+1,1); fxk(1) = c'*xk;
KKT_lk = fxk;KKT_lk(1) = norm(Ax(xk,p,q)-b);
KKT_xk = fxk;KKT_xk(1) = norm(xk-prox(xk-c-Aty(lk(1:m+n),p,q)));
%% 2. Main Loop
res0 = max(KKT_xk(1),KKT_lk(1));tic;
disp(['Let us go :) ======= ',datestr(now)]);
aa = 'A-ADMM: it = %6d, KKT(xk) = %4.2e, KKT(lk) = %4.2e, fk = %4.5e\n';
bb = '        it = %6d, KKT(xk) = %4.2e, KKT(lk) = %4.2e, fk = %4.5e\n';

tic;
for k = 1 : maxit
    resk = max(KKT_xk(k),KKT_lk(k));
    if k == 1 || k == maxit
        if k == 1
            fprintf(aa,k-1,KKT_xk(k),KKT_lk(k),fxk(k));
        else
            fprintf(bb,k-1,KKT_xk(k),KKT_lk(k),fxk(k));
        end
    else
        if res0/resk >= 2
            res0 = resk;fprintf(bb,k-1,KKT_xk(k),KKT_lk(k),fxk(k));
        end
    end
    % 2.1 Iterate
    ak = bk;bk1 = bk/(1+ak);
    gk1 = (gk+muf*ak)/(1+ak);
    etafk = (1+ak)*gk + muf*ak;
    sgk = 1/bk1;etagk = (1+ak)*bk;
    
    wwk = (ak*pik + wk)/(1+ak);
    wxk = (ak*gk*vk + (gk+muf*ak)*xk)/etafk;
    
    hlk = lk - 1/bk*[Ax(xk,p,q)-b;xk-wk] + ak/bk*[z0;-(pik-wk)];
    cAw = - Atb - wk;cAlk = Aty(hlk(1:m+n),p,q)+hlk(m+n+1:end);
    dd = etafk*wxk - ak^2*(c+cAlk+sgk*cAw);
    % 2.2 Update
    tt = sgk*ak^2;sg = 1+etafk/tt;
    xk1 = (dd-Aty(invAAt(Ax(dd,p,q),p,q,sg),p,q))/(etafk+tt);
    vk1 = xk1 + (xk1 - xk)/ak;
    blk = lk + ak/bk*[Ax(vk1,p,q)-b;vk1-pik];
    wk1 = prox(wwk-ak^2/etagk*(-blk(m+n+1:end)));
    pik1 = wk1 + (wk1-wk)/ak;
    lk1 = lk + ak/bk*[Ax(vk1,p,q)-b;vk1-pik1];
    % 2.3 Transfer
    gk = gk1;bk = bk1;xk = xk1;vk = vk1;wk = wk1;pik = pik1;lk = lk1;
    
    fxk(k+1) = c'*xk;KKT_lk(k+1) = norm(Ax(xk,p,q)-b);
    KKT_xk(k+1) = norm(xk-prox(xk-c-Aty(lk(1:m+n),p,q)));
    % 2.4 Check residual
    rr = [KKT_xk(k+1)/(1+KKT_xk(1)),KKT_lk(k+1)/(1+KKT_lk(1))];
    if max(rr) <= KKT_Tol
        gg = 'CONV at it = %3d, KKT(xk) = %4.2e, KKT(lk) = %4.2e, time used %4.2es\n';
        fprintf(1,gg,k,rr,toc);break;
    else
        if k == maxit
            gg = 'NOTCONV at it = %3d, KKT(xk) = %4.2e, KKT(lk) = %4.2e, time used %4.2es\n';
            fprintf(2,gg,k,rr,toc);
        end
    end
end
%% 3. Plot
close all;figure;
KKT_xk = KKT_xk(KKT_xk>-inf);KKT_xk = KKT_xk(2:end);KKT_xk = KKT_xk/(1+KKT_xk(1));
KKT_lk = KKT_lk(KKT_lk>-inf);KKT_lk = KKT_lk(2:end);KKT_lk = KKT_lk/(1+KKT_lk(1));
fxk = fxk(fxk>-inf);efxk = abs(fxk-fxk(end));efxk = efxk(2:end);efxk = efxk/(1+efxk(1));

loglog(1:length(KKT_xk),KKT_xk,'b-','LineWidth',1.5);hold on;
loglog(1:length(KKT_lk),KKT_lk,'k-','LineWidth',1.5);
loglog(1:length(efxk),efxk,'r-','LineWidth',1.5);hold off;

xlabel('$\log k$','Interpreter', 'latex','FontSize',14);
ylabel('$\log {\rm Errors}$','Interpreter', 'latex','FontSize',14);
set(legend('${\rm KKT}(x_k)$','${\rm KKT}(\lambda_k)$',...
    '$|f(x_k)-f^*|$'),'Interpreter','latex','FontSize',14);