%% Accelerated ADMM for problem CLASS 1:  Transport-like problem
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
addpath '../AMG';
%% 0. Problem Setting

% m = 5e2;n = m;

% prob = 0;       % half l1-norm = 0.5*norm(r-l,1)
% disp(['  0.5*||l-r||_1 = ', num2str(0.5*norm(l-r,1),'%4.5e')]);
% prob = 1;  % Assignment problem
prob = 2;  % Optimal transport
% prob = 3;  % Capacity constrained mass transport

% [~,~,~,p,q,gama] = prob_set(m,n,prob);

load('./InputData/data1-1000.mat');


%% 1. Algorithm Setting
% 1.1 auxiliary vars
prox = @(x) min(max(0,x),gama);
b = [r;l];Atb = Aty(b,p,q);z0 = zeros(m+n,1);
% 1.2 parameters
maxit = 1e2;  KKT_Tol = 1e-6; bk = 1;
SsN_IT = 50; SsN_Tol1 = 1e-11; nu = 0.2;delta = 0.9; ll_max = 500;
% 1.3 initial guess
%%%%%%%%%% random initial %%%%%%%%%%
% xk = rand(m*n,1);vk = xk;lk = rand(m+n,1);
% xk = spalloc(m*n,1,0); vk = xk;lk = zeros(m+n,1);
% xk = c; vk = xk;lk = b;
%%%%%%%%%% MATLAB linprog %%%%%%%%%%
% % options = optimoptions('linprog','Algorithm','interior-point-legacy',...
% %     'Display','iter','MaxIterations',100);
% options = optimoptions('linprog','Algorithm','interior-point-legacy',...
%     'Display','iter','MaxIterations',10,'OptimalityTolerance',1e0);
% A = [kron(speye(n),p');kron(q',speye(m))];
% [mat_xk,fk,~,~,mat_lk] = linprog(c,[],[],A,b,spalloc(m*n,1,0),gama,[],options);
% mat_lk = mat_lk.eqlin;xk = mat_xk; vk = xk; lk = mat_lk;
% % disp(['  fxk = ', num2str(fk,'%4.5e')]);
%%%%%%%%%% MATLAB linprog %%%%%%%%%%
%%%%%%%%%% A-ADMM  %%%%%%%%%%
if prob > 0
    % if inner_solver == 2
    res_0 = 0; maxit_0 = 1e2;
else
    res_0 = 5e-2;maxit_0 = inf; % This is for problem 3
end
[xk0,lk0] = warmup_class1(c,r,l,p,q,gama,res_0,maxit_0);
xk = xk0;vk = xk;lk = lk0;
%%%%%%%%%% A-ADMM  %%%%%%%%%%
% 1.4 residuals
fxk = - inf*ones(maxit+1,1); fxk(1) = c'*xk;
KKT_lk = fxk;KKT_lk(1) = norm(Ax(xk,p,q)-b);
KKT_xk = fxk;KKT_xk(1) = norm(xk-prox(xk-c-Aty(lk,p,q)));
% 1.5 inner solver %%%%%%%%%%
% inner_solver = 1; % direct solver
% inner_solver = 2; % PCG
% inner_solver = 3; % Aug_PCG
inner_solver = 4; % AMG
% inner_solver = 5; % Twogrid
%% 2. Main Loop
res0 = max(KKT_xk(1),KKT_lk(1));flag_apd = 0;
tic; disp(['Let us go :) ======= ',datestr(now)]);
aa = 'APD: it = %3d, KKT(xk) = %4.2e, KKT(lk) = %4.2e, fk = %4.5e\n';
switch inner_solver
    case 1   % direct solver
        hh = '===Direct: (%1d,%1d), its = %1d, rel_res = %1d';
    case 2   % PCG
        % Be careful with the stagnation
        pcg_options = struct('retol',1e-11,'maxit',1e4,'precd',2,'guess',[]);
        hh = '===PCG: (%1d,%1d), its = %4d, rel_res = %4.2e';
    case 3 % Be careful with the stagnation
        pcg_options = struct('retol',1e-11,'maxit',1e4,'precd',2,'guess',[]);
        hh = '===Aug_PCG: (%3d,%1d), its = %4d, rel_res = %4.2e';
    case {4,5}   % AMG or Twogrid
        amg_options = struct('retol',1e-11,'bigph',1,'maxit',30,'theta',1/4,...
            'smoth',5,'cycle','w','isnsp',1,'inter',1,'guess',[]);
        hh = '===Hybrid_AMG: (%3d,%1d), its = %4d, rel_res = %4.2e';
end
cc = ['   SsN: it = %3d, |Fk| = %4.2e, ll = %3d',hh,'\n'];
bb = ['   SsN: it = %3d, |Fk| = %4.2e, ll = %3d',hh,'\n\n'];

SumAMG=0;
TotalAMG=0;
FailAMG=0;
MaxAMG=0;


tic;SsN_itnum = nan(maxit,1);PCG_itnum = nan(maxit,3);
for k = 1 : maxit
    resk = max(KKT_xk(k),KKT_lk(k));
    if k == 1 || k == maxit
        flag_apd = 1; fprintf(aa,k-1,KKT_xk(k),KKT_lk(k),fxk(k));
    else
        if res0/resk >= 0
            res0 = resk;flag_apd = 1;
            fprintf(aa,k-1,KKT_xk(k),KKT_lk(k),fxk(k));
        end
    end
    % 2.1 iterate
    %                             ak = 1+rand;
    ak = sqrt(k^(2)*bk);
%     ak=1;
%     if k<=10
%         ak=1;
%     else
%         ak=0.4;
%     end
    bk1 = bk/(1+ak);tk  = bk*(1+ak)/ak^2;
    
%     SsN_Tol=SsN_Tol1;
    SsN_Tol=max(bk1/(k^2),SsN_Tol1);
    
    wk  = -c + bk*(xk+ak*vk)/ak^2;
    wlk = bk1*(lk-1/bk*(Ax(xk,p,q)-b))-b;
    % 2.2 SsN
    ssn_it = 0;lk_new = lk;
    zk = 1/tk*(wk-Aty(lk_new,p,q));
    Fk_new = bk1*lk_new - Ax(prox(zk),p,q) - wlk;
    Fk_res = norm(Fk_new);
    
    pcg_itnum = -inf(SsN_IT,1);flag_ssn = flag_apd;
    if flag_ssn
        fprintf('   SsN: it = %3d, |Fk| = %4.2e\n',0,Fk_res);
    end
    while norm(Fk_new) > SsN_Tol
        ssn_it = ssn_it + 1;lk_old = lk_new;
        zk = 1/tk*(wk-Aty(lk_old,p,q));
        s = (zk>=0) & (zk<=gama); t = zeros(m+n,1);
        % Jk*zeta = -Fk(lk_old)
        T = spdiags(t,0,length(t),length(t));H0 = ASAt(s,p,q);
        % Jk = bk1*speye(m+n) + 1/tk*(diag(t)+A*diag(s)*A');
        Fk_old = bk1*lk_old - Ax(prox(zk),p,q) - wlk;
        switch inner_solver
            case 1
                Jk = bk1*speye(m+n)+(T+H0)/tk;
                zeta = Jk \ (-Fk_old);info = [0,0];itpcg = 1;respcg = 0;
            case 2
                if pcg_options.precd == 5; pcg_options.nf = n; end
                Jk = bk1*speye(m+n)+(T+H0)/tk;
                info = [0,0];[zeta,itpcg,respcg] = PCG(Jk,-Fk_old,pcg_options);
            case {3,4,5}
                prob_data.bk1 = bk1;prob_data.tk = tk;
                prob_data.q = q;prob_data.p = p;prob_data.T = T;
                prob_data.H0 = H0;prob_data.z = - Fk_old;
                if inner_solver == 3
                    [zeta,itpcg,respcg,info] = aug_PCG(prob_data,pcg_options);
                end
                if inner_solver == 4
                    [zeta,itpcg,respcg,info] = Hybrid_AMG(prob_data,amg_options);
                    
                    if itpcg == amg_options.maxit
%                         warning('AMG failed!!! respcg = %f',respcg);
                        FailAMG=FailAMG+1;
%                         save('wrongAMG.mat','prob_data','amg_options');
                    else
                        MaxAMG=max(MaxAMG,itpcg);
                    end
                    if itpcg>0
                        TotalAMG=TotalAMG+1;
                    end
                    
                    
                end
                if inner_solver == 5
                    [zeta,itpcg,respcg,info] = Hybrid_twogrid(prob_data,amg_options);
                end
        end
        pcg_itnum(ssn_it) = itpcg;
        % 2.3 Line search
        f0 = bk1/2*norm(lk_old)^2 - wlk'*lk_old;
        if prob < 3
            cFk_old = f0 + 0.5*tk*norm(prox(zk))^2;
        else
            cFk_old = f0 + 0.5*tk*(norm(zk)^2 - norm(zk-prox(zk))^2);
        end
        
        ll = 0;lk_new = lk_old + delta^ll*zeta;
        f0 = bk1/2*norm(lk_new)^2 - wlk'*lk_new;
        zk = 1/tk*(wk-Aty(lk_new,p,q));
        if prob < 3
            cFk_new = f0 + 0.5*tk*norm(prox(zk))^2;
        else
            cFk_new = f0 + 0.5*tk*(norm(zk)^2 - norm(zk-prox(zk))^2);
        end
        
        ress = abs(Fk_old'*zeta);
        while cFk_new > cFk_old -  nu*delta^ll * ress
            ll = ll + 1;lk_new = lk_old + delta^ll*zeta;
            f0 = bk1/2*norm(lk_new)^2 - wlk'*lk_new;
            zk = 1/tk*(wk-Aty(lk_new,p,q));
            if prob < 3
                cFk_new = f0 + 0.5*tk*norm(prox(zk))^2;
            else
                cFk_new = f0 + 0.5*tk*(norm(zk)^2 - norm(zk-prox(zk))^2);
            end
            if ll == ll_max
                break;
            end
        end
        Fk_new = bk1*lk_new - Ax(prox(zk),p,q) - wlk;
        if norm(Fk_new) <= SsN_Tol
            if flag_ssn
                fprintf(bb,ssn_it,norm(Fk_new),ll,info,itpcg,respcg);
            end
            break;
        else
            if abs(norm(Fk_old)-norm(Fk_new)) < SsN_Tol/100
                if flag_ssn
                    fprintf(bb,ssn_it,norm(Fk_new),ll,info,itpcg,respcg);
                end
                break;
            end
        end
        if ssn_it == SsN_IT
            if flag_ssn
                fprintf(bb,ssn_it,norm(Fk_new),ll,info,itpcg,respcg);
            end
            break;
        end
        if Fk_res/norm(Fk_new) >= 2
            Fk_res = norm(Fk_new);
            if flag_ssn
                fprintf(cc,ssn_it,norm(Fk_new),ll,info,itpcg,respcg);
            end
        end
    end
    flag_apd = 0;lk1 = lk_new;xk1 = prox(zk);vk1 = xk1 + (xk1-xk)/ak;
    %%%Restart%%%%%%%%%%%%%%%%%%%%%%%%%
    KKT_lk(k+1) = norm(Ax(xk1,p,q)-b);
    KKT_xk(k+1) = norm(xk1-prox(xk1-c-Aty(lk1,p,q)));
    
    rr = [KKT_xk(k+1)/(1+KKT_xk(1)),KKT_lk(k+1)/(1+KKT_lk(1))];
    if bk1 < 1e-8 && max(rr) > resk
        xk1 = xk;lk1 = lk;vk1 = xk;bk1 = rand;
        fprintf(2,['==============Restart at it = ',...
            num2str(k),'==============','\n']);
    end
    %%%Restart%%%%%%%%%%%%%%%%%%%%%%%%%
    bk = bk1;xk = xk1;lk = lk1;vk = vk1;
    
    fxk(k+1) = c'*xk;KKT_lk(k+1) = norm(Ax(xk,p,q)-b);
    KKT_xk(k+1) = norm(xk-prox(xk-c-Aty(lk(1:m+n),p,q)));
    
    SsN_itnum(k) = ssn_it;
    if ssn_it
        pcg_itnum = pcg_itnum(1:ssn_it);re = fix(sum(pcg_itnum)/ssn_it);
        PCG_itnum(k,:) = [min(pcg_itnum),re,max(pcg_itnum)];
        
        SumAMG=SumAMG+sum(pcg_itnum);
        
    end
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
% figure(1); clf;
KKT_xk = KKT_xk(KKT_xk>-inf);KKT_xk = KKT_xk(2:end);KKT_xk = KKT_xk/(1+KKT_xk(1));
KKT_lk = KKT_lk(KKT_lk>-inf);KKT_lk = KKT_lk(2:end);KKT_lk = KKT_lk/(1+KKT_lk(1));
fxk = fxk(fxk>-inf);efxk = abs(fxk-fxk(end));efxk = efxk(3:end);efxk = efxk/(1+efxk(1));

SsN_itnum = SsN_itnum(1:k);PCG_itnum = PCG_itnum(1:k,:);
figure(1)
subplot(1,3,1);
loglog(1:length(efxk),efxk,'dk-','LineWidth',1.5);hold on;
loglog(1:length(KKT_xk),KKT_xk,'^b-','LineWidth',1.5);
loglog(1:length(KKT_lk),KKT_lk,'or-','LineWidth',1.5); hold off;

xlabel('$k$','Interpreter', 'latex','FontSize',14);
ylabel('$\log {\rm Errors}$','Interpreter', 'latex','FontSize',14);
set(legend('$|f(x_k)-f^*|$','${\rm KKT}(x_k)$',...
    '${\rm KKT}(\lambda_k)$'),'Interpreter','latex','FontSize',13);

subplot(1,3,2);
semilogy(1:length(SsN_itnum),SsN_itnum','^b-','LineWidth',1.5);

xlabel('$k$','Interpreter', 'latex','FontSize',14);
ylabel('$\#{\rm SsN}$','Interpreter', 'latex','FontSize',14);
set(legend('${\rm SsN}(k)$'),'Interpreter','latex','FontSize',14);

subplot(1,3,3);
semilogy(1:length(PCG_itnum(:,3)),PCG_itnum(:,3)','or-','LineWidth',1.5);hold on;
semilogy(1:length(PCG_itnum(:,2)),PCG_itnum(:,2)','^b-','LineWidth',1.5);
semilogy(1:length(PCG_itnum(:,1)),PCG_itnum(:,1)','dk-','LineWidth',1.5);hold off;

xlabel('$k$','Interpreter', 'latex','FontSize',14);
ylabel('$\#{\rm AMG}$','Interpreter', 'latex','FontSize',14);
set(legend('${\rm max\_AMG}(k)$','${\rm ava\_AMG}(k)$','${\rm min\_AMG}(k)$'),...
    'Interpreter','latex','FontSize',14);

figure(2)
semilogy(1:length(PCG_itnum(:,3)),PCG_itnum(:,3)','or-','LineWidth',1.5);hold on;
semilogy(1:length(PCG_itnum(:,2)),PCG_itnum(:,2)','^b-','LineWidth',1.5);
semilogy(1:length(PCG_itnum(:,1)),PCG_itnum(:,1)','dk-','LineWidth',1.5);hold off;

xlabel('$k$','Interpreter', 'latex','FontSize',14);
ylabel('$\#{\rm iter}$','Interpreter', 'latex','FontSize',14);
set(legend('${\rm max\_AMG}(k)$','${\rm ava\_AMG}(k)$','${\rm min\_AMG}(k)$'),...
    'Interpreter','latex','FontSize',14);


figure(3)
loglog(1:length(KKT_xk),KKT_xk,'^b-','LineWidth',1.5); hold on;
loglog(1:length(KKT_lk),KKT_lk,'or-','LineWidth',1.5); hold off;

xlabel('$k$','Interpreter', 'latex','FontSize',14);
ylabel('$\log {\rm Errors}$','Interpreter', 'latex','FontSize',14);
set(legend('${\rm KKT}(x_k)$',...
    '${\rm KKT}(\lambda_k)$'),'Interpreter','latex','FontSize',13);




