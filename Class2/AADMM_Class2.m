%% Accelerated ADMM for problem CLASS 2: Partial OT
%        min_{x>=0,y>=0,z>=0}  <c,x>
%             s.t. G*x + IY*y + IZ*z = b
% where c in R^{mn}_+, r in R^n_+, l in R^m_+ and
%      [ In  ox 1m']        [ In]        [ O ]        [ r ]
%  G = [ 1n' ox Im ]   IY = [ O ]   IZ = [ Im]    b = [ l ]
%      [   1mn'    ]        [ O ]        [ O ]        [ mu]
% with mu in (0, min(<r,1n>, <l,1m>)].

addpath '..';
%% 0. Problem Setting
% m = 40e2;n = m;
% C = rand(m,n); c = C(:);
% r = rand(n,1);q = ones(n,1);
% l = rand(m,1);p = ones(m,1);
% mu = rand*min(r'*q,l'*p);
% phi = ones(m*n,1);

load('./InputData/data4-1000.mat');








%% 1. Algorithm Setting
% 1.1 auxiliary vars
prox = @(x) max(0,x);b = [r;l;mu];wc = [c;zeros(n+m,1)];
Htb = [Aty(b(1:m+n),p,q)+b(end)*phi;b(1:n+m)];z0 = zeros(m+n+1,1);
% 1.2 parameters
muf = 0; gk = 1;bk = 1; maxit = 5e3; KKT_Tol = 1e-6;
% 1.3 initial guess
%%%%%%%%%% MATLAB linprog %%%%%%%%%%
% % options = optimoptions('linprog','Algorithm','interior-point-legacy',...
% %     'Display','iter','MaxIterations',100,'OptimalityTolerance',1e-10);
% options = optimoptions('linprog','Algorithm','interior-point-legacy',...
%     'Display','iter','MaxIterations',10,'OptimalityTolerance',1e-2);
% A = [kron(speye(n),p');kron(q',speye(m))];G = [A;phi'];
% IY = [speye(n);spalloc(m,n,0);spalloc(1,n,0)];
% IZ = [spalloc(n,m,0);speye(m);spalloc(1,m,0)]; 
% H = [G,IY,IZ]; x0 = H(1,:)';
% [mat_uk,fk,~,~,mat_lk] = linprog(wc,[],[],H,b,x0,[],[],options);
% mat_lk = mat_lk.eqlin;uk = mat_uk; vk = uk; wk = uk;pik = wk;lk = mat_lk; 
% disp(['  fxk = ', num2str(fk,'%4.5e')]); 
%%%%%%%%%% MATLAB linprog %%%%%%%%%%
%%%%%%%%%% random initial %%%%%%%%%%
uk = spalloc(m*n+n+m,1,0); vk = uk;wk = uk;pik = wk;lk = [z0;uk];
% uk = [c;r;l]; vk = uk;wk = uk;pik = wk;lk = [b;uk];
% 1.4 residuals
xk = uk(1:m*n);yk = uk(m*n+1:m*n+n);zk = uk(m*n+n+1:end);
fxk = - inf*ones(maxit+1,1); fxk(1) = c'*xk;
KKT_lk = fxk;KKT_lk(1) = norm([Ax(xk,p,q)+[yk;zk];phi'*xk]-b);
KKT_zk = fxk;KKT_zk(1) = norm(zk-max(zk-lk(n+1:n+m),0));
KKT_yk = fxk;KKT_yk(1) = norm(yk-max(yk-lk(1:n),0));
KKT_xk = fxk;KKT_xk(1) = norm(xk-max(xk-c-(Aty(lk(1:m+n),p,q)+lk(m+n+1)*phi),0));
%% 2. Main Loop
res0 = max(max(max(KKT_xk(1),KKT_yk(1)),KKT_zk(1)),KKT_lk(1));
tic; disp(['Let us go :) ======= ',datestr(now)]);
aa = 'KKT(xk) = %4.2e, KKT(yk) = %4.2e, KKT(zk) = %4.2e, KKT(lk) = %4.2e';
dd1 = ['A-ADMM: it = %6d, ',aa,', fk = %4.5e\n'];
dd2 = ['        it = %6d, ',aa,', fk = %4.5e\n'];

tic;
for k = 1 : maxit
    resk = max(max(max(KKT_xk(k),KKT_yk(k)),KKT_zk(k)),KKT_lk(k));
    if k == 1 || k == maxit
        if k == 1
            fprintf(dd1,k-1,KKT_xk(k),KKT_yk(k),KKT_zk(k),KKT_lk(k),fxk(k));
        else
            fprintf(dd2,k-1,KKT_xk(k),KKT_yk(k),KKT_zk(k),KKT_lk(k),fxk(k));
        end
    else
        if res0/resk >= 2
            res0 = resk;
            fprintf(dd2,k-1,KKT_xk(k),KKT_yk(k),KKT_zk(k),KKT_lk(k),fxk(k));
        end
    end
    % 1.1 Iterate
    ak = bk;bk1 = bk/(1+ak);
    gk1 = (gk+muf*ak)/(1+ak);
    etafk = (1+ak)*gk + muf*ak;
    sgk = 1/bk1;etagk = (1+ak)*bk;
    
    wwk = (ak*pik + wk)/(1+ak);
    wuk = (ak*gk*vk + (gk+muf*ak)*uk)/etafk;
    
    hlk = lk - 1/bk*[[Ax(xk,p,q)+[yk;zk];phi'*xk]-b;uk-wk] + ak/bk*[z0;-(pik-wk)];
    cAw = - Htb - wk;cAlk = hlk(m+n+2:end);
    cAlk = cAlk + [Aty(hlk(1:m+n),p,q)+hlk(m+n+1)*phi;hlk(1:n+m)];
    dd = etafk*wuk - ak^2*(wc+cAlk+sgk*cAw);
    % 1.2 Update
    tt = sgk*ak^2;sg = 1+etafk/tt;
    Hdd = [Ax(dd(1:m*n),p,q) + dd(m*n+1:end);phi'*dd(1:m*n)];
    ff = invHHt(Hdd,p,q,sg,phi);
    uk1 = (dd-[Aty(ff(1:m+n),p,q)+ff(end)*phi;ff(1:m+n)])/(etafk+tt);
    
    vk1 = uk1 + (uk1 - uk)/ak;
    b0 = [Ax(vk1(1:m*n),p,q)+vk1(m*n+1:end);phi'*vk1(1:m*n)]-b;
    blk = lk + ak/bk*[b0;vk1-pik];
    wk1 = prox(wwk-ak^2/etagk*(-blk(m+n+2:end)));
    
    pik1 = wk1 + (wk1-wk)/ak;
    lk1 = lk + ak/bk*[b0;vk1-pik1];
    % 1.3 Transfer
    gk = gk1;bk = bk1;uk = uk1;vk = vk1;
    wk = wk1;pik = pik1;lk = lk1;
    
    xk = uk(1:m*n);yk = uk(m*n+1:m*n+n);
    zk = uk(m*n+n+1:end);fxk(k+1) = c'*xk;
    KKT_lk(k+1) = norm([Ax(xk,p,q)+[yk;zk];phi'*xk]-b);
    KKT_zk(k+1) = norm(zk-max(zk-lk(n+1:n+m),0));
    KKT_yk(k+1) = norm(yk-max(yk-lk(1:n),0));
    KKT_xk(k+1) = norm(xk-max(xk-c-(Aty(lk(1:m+n),p,q)+lk(m+n+1)*phi),0));
    % 2.4 Check residual          
    rr1 = [KKT_xk(k+1)/(1+KKT_xk(1)),KKT_yk(k+1)/(1+KKT_yk(1))];
    rr = [rr1 KKT_zk(k+1)/(1+KKT_zk(1)),KKT_lk(k+1)/(1+KKT_lk(1))];
    if max(rr) <= KKT_Tol
        fprintf(1,['CONV at it = %3d, ',aa,', time used %4.2es\n'],k,rr,toc);break;
    else
        if k == maxit
            fprintf(2,['NOTCONV at it = %3d, ',aa,', time used %4.2es\n'],k,rr,toc);
        end
    end
end
%%
close all;figure;
KKT_zk = KKT_zk(KKT_zk>-inf);KKT_zk = KKT_zk(2:end);KKT_zk = KKT_zk/(1+KKT_zk(1));
KKT_yk = KKT_yk(KKT_yk>-inf);KKT_yk = KKT_yk(2:end);KKT_yk = KKT_yk/(1+KKT_yk(1));
KKT_xk = KKT_xk(KKT_xk>-inf);KKT_xk = KKT_xk(2:end);KKT_xk = KKT_xk/(1+KKT_xk(1));
KKT_lk = KKT_lk(KKT_lk>-inf);KKT_lk = KKT_lk(2:end);KKT_lk = KKT_lk/(1+KKT_lk(1));
fxk = fxk(fxk>-inf);efxk = abs(fxk-fxk(end));efxk = efxk(2:end);efxk = efxk/(1+efxk(1));

loglog(1:length(KKT_lk),KKT_lk,'k-','LineWidth',1.5);hold on;
loglog(1:length(KKT_xk),KKT_xk,'m-','LineWidth',1.5);
loglog(1:length(KKT_yk),KKT_yk,'g-','LineWidth',1.5);
loglog(1:length(KKT_zk),KKT_zk,'b-','LineWidth',1.5);
loglog(1:length(efxk),efxk,'r-','LineWidth',1.5);hold off;

xlabel('$k$','Interpreter', 'latex','FontSize',14);
ylabel('$\log {\rm Errors}$','Interpreter', 'latex','FontSize',14);
set(legend('${\rm KKT}(\lambda_k)$','${\rm KKT}(x_k)$','${\rm KKT}(y_k)$',...
    '${\rm KKT}(z_k)$','$|f(x_k)-f^*|$'),'Interpreter','latex','FontSize',14);