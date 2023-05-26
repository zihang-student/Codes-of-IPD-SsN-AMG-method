%% Inexact SsN-AMG-APD method for problem CLASS 2: Partial OT
%        min_{x>=0,y>=0,z>=0}  <c,x>
%             s.t. G*x + IY*y + IZ*z = b
% where c in R^{mn}_+, r in R^n_+, l in R^m_+ and
%      [ In  ox 1m']        [ In]        [ O ]        [ r ]
%  G = [ 1n' ox Im ]   IY = [ O ]   IZ = [ Im]    b = [ l ]
%      [   1mn'    ]        [ O ]        [ O ]        [ mu]
% with mu in (0, min(<r,1n>, <l,1m>)].

addpath '..'; 
addpath '../AMG';
%% 0. Problem Setting
% m = 5e2;n = m;
% C = rand(m,n); c = C(:);
% r = rand(n,1);q = ones(n,1);
% l = rand(m,1);p = ones(m,1);
% mu = rand*min(r'*q,l'*p);
% phi = ones(m*n,1);

load('./InputData/data4-1000.mat');


%% 1. Algorithm Setting
% 1.1 auxiliary vars
prox = @(x) max(0,x);b = [r;l;mu];wc = [c;zeros(n+m,1)];
% 1.2 parameters
maxit = 1e2; KKT_Tol = 1e-6; bk = 1;
SsN_IT = 50; SsN_Tol1 = 1e-10; nu = 0.2;delta = 0.9; ll_max = 500;
% 1.3 initial guess
%%%%%%%%%% random initial %%%%%%%%%%
% uk = rand(m*n+n+m,1); vk = uk; lk = rand(m+n+1,1);
% uk = spalloc(m*n+n+m,1,0); vk = uk; lk = zeros(m+n+1,1);
% uk = [c;r;l]; vk = uk;lk = b;
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
% mat_lk = mat_lk.eqlin;uk = mat_uk; vk = uk; lk = mat_lk;
% disp(['  fxk = ', num2str(fk,'%4.5e')]);
%%%%%%%%%% MATLAB linprog %%%%%%%%%%
%%%%%%%%%% A-ADMM  %%%%%%%%%%
res_0 = 0;maxit_0 = 1e2;
% res_0 = 1e-1;maxit_0 = inf;
[uk0,lk0] = warmup_class2(c,r,l,p,q,mu,phi,res_0,maxit_0);
uk = uk0;lk = lk0; vk = uk;
%%%%%%%%%% A-ADMM  %%%%%%%%%%
% 1.4 residuals
xk = uk(1:m*n);yk = uk(m*n+1:m*n+n);zk = uk(m*n+n+1:end);
fxk = - inf*ones(maxit+1,1); fxk(1) = c'*xk;
KKT_lk = fxk;KKT_lk(1) = norm([Ax(xk,p,q)+[yk;zk];phi'*xk]-b);
KKT_zk = fxk;KKT_zk(1) = norm(zk-max(zk-lk(n+1:n+m),0));
KKT_yk = fxk;KKT_yk(1) = norm(yk-max(yk-lk(1:n),0));
KKT_xk = fxk;KKT_xk(1) = norm(xk-max(xk-c-(Aty(lk(1:m+n),p,q)+lk(m+n+1)*phi),0));
% 1.5 inner solver %%%%%%%%%%
% inner_solver = 1; % direct solver
% inner_solver = 2; % PCG
% inner_solver = 3; % Aug_PCG
inner_solver = 4; % AMG
% inner_solver = 5; % Twogrid
%% 2. Main Loop
res0 = max(max(max(KKT_xk(1),KKT_yk(1)),KKT_zk(1)),KKT_lk(1));flag_apd = 0;
tic; disp(['Inexact SsN-AMG-APD: Let us go :) ======= ',datestr(now)]);
aa = 'KKT(xk) = %4.2e, KKT(yk) = %4.2e, KKT(zk) = %4.2e, KKT(lk) = %4.2e';
switch inner_solver
    case 1   % direct solver
        hh = '===Direct: (%1d,%1d), its = %1d, rel_res = %1d';
    case 2   % PCG
        pcg_options = struct('retol',1e-11,'maxit',1e4,'precd',2,'guess',[]);
        hh = '===PCG: (%1d,%1d), its = %4d, rel_res = %4.2e';
    case 3   % Aug_PCG
        pcg_options = struct('retol',1e-11,'maxit',1e4,'precd',2,'guess',[]);
        hh = '===Aug_PCG: (%3d,%1d), its = %4d, rel_res = %4.2e';
    case {4,5}   % AMG or Twogrid
        amg_options = struct('retol',1e-11,'bigph',1,'maxit',4e1,'theta',1/4,...
            'smoth',10,'cycle','w','isnsp',1,'inter',1,'guess',[]);
        hh = '===Hybrid_AMG: (%3d,%1d), its = %4d, rel_res = %4.2e';
end
cc = ['   SsN: it = %3d, |Fk| = %4.2e, ll = %3d',hh,'\n'];
bb = ['   SsN: it = %3d, |Fk| = %4.2e, ll = %3d',hh,'\n\n'];


SumAMG=0;
TotalAMG=0;
FailAMG=0;
MaxAMG=0;


tic;SsN_itnum = -inf(maxit,1);PCG_itnum = -inf(maxit,3);
for k = 1 : maxit
    resk = max(max(max(KKT_xk(k),KKT_yk(k)),KKT_zk(k)),KKT_lk(k));
    if k == 1 || k == maxit
        flag_apd = 1;dd = ['APD: it = %3d, ',aa,', fk = %4.5e\n'];
        fprintf(dd,k-1,KKT_xk(k),KKT_yk(k),KKT_zk(k),KKT_lk(k),fxk(k));
    else
        if res0/resk >= 2
            res0 = resk;flag_apd = 1;
            dd = ['APD: it = %3d, ',aa,', fk = %4.5e\n'];
            fprintf(dd,k-1,KKT_xk(k),KKT_yk(k),KKT_zk(k),KKT_lk(k),fxk(k));
        end
    end
    % 2.1 iterate
%         ak = 1+rand;

%     ak=0.4;
%     if k<=7
%         ak=1;
%     else
%         ak=0.4;
%     end
    ak = sqrt(k^2*bk);
    bk1 = bk/(1+ak);tk = bk*(1+ak)/ak^2;
    
    SsN_Tol=max(bk1/(k^2),SsN_Tol1);
    
    wk  = -wc + bk*(uk+ak*vk)/ak^2;
    wlk = bk1*(lk-1/bk*([Ax(xk,p,q)+[yk;zk];phi'*xk]-b))-b;
    % 2.2 SsN
    ssn_it = 0;lk_new = lk;
    Htlk = [Aty(lk_new(1:m+n),p,q)+lk_new(m+n+1)*phi;lk_new(1:m+n)];
    
    zk = 1/tk*(wk - Htlk);pzk = prox(zk);
    ppk = Ax(pzk(1:m*n),p,q)+[pzk(m*n+1:m*n+n);pzk(m*n+n+1:end)];
    Hpzk = [ppk;phi'*pzk(1:m*n)];
    Fk_new = bk1*lk_new - Hpzk - wlk;Fk_res = norm(Fk_new);
    
    pcg_itnum = -inf(SsN_IT,1);flag_ssn = flag_apd;
    if flag_ssn
        fprintf('   SsN: it = %3d, |Fk| = %4.2e\n',0,Fk_res);
    end
    while norm(Fk_new) > SsN_Tol
        ssn_it = ssn_it + 1;lk_old = lk_new;
        Htlk = [Aty(lk_old(1:m+n),p,q)+lk_old(m+n+1)*phi;lk_old(1:m+n)];
        zk = 1/tk*(wk - Htlk); s = zk(1:m*n)>=0;t = zk(m*n+1:end)>=0;
        % Jk*zeta = -Fk(lk_old)
        % Jk = bk1*speye(m+n+1) + 1/tk*(cT + cH0);
        % cT = [diag(t) zeros(m+n,1);zeros(1,m+n+1)]
        %                      [A*diag(s)*A'       A*diag(s)*phi ]
        % cH0 = G*diag(s)*G' = [                                 ]
        %                      [phi'*diag(s)*A'  phi'*diag(s)*phi]
        T = spdiags(t,0,length(t),length(t));H0 = ASAt(s,p,q);
        
        pzk = prox(zk);
        ppk = Ax(pzk(1:m*n),p,q)+[pzk(m*n+1:m*n+n);pzk(m*n+n+1:end)];
        Hpzk = [ppk;phi'*pzk(1:m*n)];Fk_old = bk1*lk_old - Hpzk - wlk;
        switch inner_solver
            case 1
                cT  = [T zeros(m+n,1);zeros(1,m+n+1)];
                ss = Ax(s.*phi,p,q);cH0 = [H0 ss; ss' phi'*(s.*phi)];
                zeta = (bk1*speye(m+n+1)+(cT+cH0)/tk) \ (-Fk_old);
                info = [0,0];itpcg = 1;respcg = 0;
            case 2
                cT  = [T zeros(m+n,1);zeros(1,m+n+1)];
                ss = Ax(s.*phi,p,q);cH0 = [H0 ss; ss' phi'*(s.*phi)];
                Jk = bk1*speye(m+n+1)+(cT+cH0)/tk;info = [0,0];
                [zeta,itpcg,respcg] = PCG(Jk,-Fk_old,pcg_options);
            case {3,4,5}
                prob_data.bk1 = bk1;prob_data.tk = tk;
                prob_data.q = q;prob_data.p = p;prob_data.s = s;
                prob_data.T = T;prob_data.H0 = H0;
                prob_data.z = - Fk_old;prob_data.phi = phi;
                if inner_solver == 3
                    [zeta,itpcg,respcg,info] = PCG4POT(prob_data,pcg_options);
                end
                if inner_solver == 4
                    [zeta,itpcg,respcg,info] = AMG4POT(prob_data,amg_options,'amg');
                    if itpcg == amg_options.maxit
                        FailAMG=FailAMG+1;
                    else
                        MaxAMG=max(MaxAMG,itpcg);
                    end
                    if itpcg>0
                        TotalAMG=TotalAMG+1;
                    end
                end
                if inner_solver == 5
                    [zeta,itpcg,respcg,info] = AMG4POT(prob_data,amg_options,'twogrid');
                    if itpcg == amg_options.maxit
                        FailAMG=FailAMG+1;
                    else
                        MaxAMG=max(MaxAMG,itpcg);
                    end
                    if itpcg>0
                        TotalAMG=TotalAMG+1;
                    end
                end
        end
        pcg_itnum(ssn_it) = itpcg;
        % 2.3 Line search
        f0 = bk1/2*norm(lk_old)^2 - wlk'*lk_old;
        cFk_old = f0 + 0.5*tk*norm(prox(zk))^2;
        
        ll = 0;lk_new = lk_old + delta^ll*zeta;
        f0 = bk1/2*norm(lk_new)^2 - wlk'*lk_new;
        Htlk = [Aty(lk_new(1:m+n),p,q)+lk_new(m+n+1)*phi;lk_new(1:m+n)];
        zk = 1/tk*(wk - Htlk);cFk_new = f0 + 0.5*tk*norm(prox(zk))^2;
        
        ress = abs(Fk_old'*zeta);
        while cFk_new > cFk_old -  nu*delta^ll * ress
            ll = ll + 1;lk_new = lk_old + delta^ll*zeta;
            f0 = bk1/2*norm(lk_new)^2 - wlk'*lk_new;
            Htlk = [Aty(lk_new(1:m+n),p,q)+lk_new(m+n+1)*phi;lk_new(1:m+n)];
            zk = 1/tk*(wk - Htlk);cFk_new = f0 + 0.5*tk*norm(prox(zk))^2;
            if ll == ll_max
                break;
            end
        end
        pzk = prox(zk);
        ppk = Ax(pzk(1:m*n),p,q)+[pzk(m*n+1:m*n+n);pzk(m*n+n+1:end)];
        Hpzk = [ppk;phi'*pzk(1:m*n)];
        Fk_new = bk1*lk_new - Hpzk - wlk;
        if norm(Fk_new) <= SsN_Tol
            if flag_ssn
                fprintf(bb,ssn_it,norm(Fk_new),ll,info,itpcg,respcg);
            end
            break;
        else
            if abs(norm(Fk_old)-norm(Fk_new)) < SsN_Tol
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
    flag_apd = 0;lk1 = lk_new;uk1 = prox(zk);vk1 = uk1 + (uk1-uk)/ak;
    %%%Restart%%%%%%%%%%%%%%%%%%%%%%%%%
    xk1 = uk1(1:m*n);yk1 = uk1(m*n+1:m*n+n);zk1 = uk1(m*n+n+1:end);
    KKT_lk(k+1) = norm([Ax(xk1,p,q)+[yk1;zk1];phi'*xk1]-b);
    KKT_zk(k+1) = norm(zk1-max(zk1-lk1(n+1:n+m),0));
    KKT_yk(k+1) = norm(yk1-max(yk1-lk1(1:n),0));
    KKT_xk(k+1) = norm(xk1-max(xk1-c-(Aty(lk1(1:m+n),p,q)+lk1(m+n+1)*phi),0));
    
    rr1 = [KKT_xk(k+1)/(1+KKT_xk(1)),KKT_yk(k+1)/(1+KKT_yk(1))];
    rr = [rr1 KKT_zk(k+1)/(1+KKT_zk(1)),KKT_lk(k+1)/(1+KKT_lk(1))];
    if bk1 < 1e-8 && max(rr) > resk
        uk1 = uk;lk1 = lk;vk1 = uk;bk1 = 10*bk1;
        fprintf(2,['==============Restart at it = ',...
            num2str(k),'==============','\n']);
    end
    %%%Restart%%%%%%%%%%%%%%%%%%%%%%%%%
    bk = bk1;uk = uk1;lk = lk1;vk = vk1;
    
    xk = uk(1:m*n);yk = uk(m*n+1:m*n+n);
    zk = uk(m*n+n+1:end);fxk(k+1) = c'*xk;
    KKT_lk(k+1) = norm([Ax(xk,p,q)+[yk;zk];phi'*xk]-b);
    KKT_zk(k+1) = norm(zk-max(zk-lk(n+1:n+m),0));
    KKT_yk(k+1) = norm(yk-max(yk-lk(1:n),0));
    KKT_xk(k+1) = norm(xk-max(xk-c-(Aty(lk(1:m+n),p,q)+lk(m+n+1)*phi),0));
    
    
    SsN_itnum(k) = ssn_it;pcg_itnum = pcg_itnum(1:ssn_it);
    PCG_itnum(k,:) = [min(pcg_itnum),fix(sum(pcg_itnum)/ssn_it),max(pcg_itnum)];
    
    SumAMG=SumAMG+sum(pcg_itnum);
    
    % 2.4 Check residual
    rr1 = [KKT_xk(k+1)/(1+KKT_xk(1)),KKT_yk(k+1)/(1+KKT_yk(1))];
    rr = [rr1 KKT_zk(k+1)/(1+KKT_zk(1)),KKT_lk(k+1)/(1+KKT_lk(1))];
    if max(rr) <= KKT_Tol
        fprintf(1,['CONV at it = %3d, ',aa,', time used %4.2es\n'],k,rr,toc);
        break;
    else
        if k == maxit
            fprintf(2,['NOTCONV at it = %3d, ',aa,', time used %4.2es\n'],k,rr,toc);
        end
    end
end
%% 3. Plot
close all;figure;
KKT_zk = KKT_zk(KKT_zk>-inf);KKT_zk = KKT_zk(2:end);KKT_zk = KKT_zk/(1+KKT_zk(1));
KKT_yk = KKT_yk(KKT_yk>-inf);KKT_yk = KKT_yk(2:end);KKT_yk = KKT_yk/(1+KKT_yk(1));
KKT_xk = KKT_xk(KKT_xk>-inf);KKT_xk = KKT_xk(2:end);KKT_xk = KKT_xk/(1+KKT_xk(1));
KKT_lk = KKT_lk(KKT_lk>-inf);KKT_lk = KKT_lk(2:end);KKT_lk = KKT_lk/(1+KKT_lk(1));
fxk = fxk(fxk>-inf);efxk = abs(fxk-fxk(end));efxk = efxk(2:end);efxk = efxk/(1+efxk(1));

SsN_itnum = SsN_itnum(1:k);PCG_itnum = PCG_itnum(1:k,:);

figure(1)
subplot(1,3,1);
semilogy(1:length(efxk),efxk,'dk-','LineWidth',1.5);hold on;
semilogy(1:length(KKT_xk),KKT_xk,'og-','LineWidth',1.5);
semilogy(1:length(KKT_yk),KKT_yk,'^m-','LineWidth',1.5);
semilogy(1:length(KKT_zk),KKT_zk,'hb-','LineWidth',1.5);
semilogy(1:length(KKT_lk),KKT_lk,'pr-','LineWidth',1.5);hold off;

xlabel('$k$','Interpreter', 'latex','FontSize',14);
ylabel('$\log {\rm Errors}$','Interpreter', 'latex','FontSize',14);
set(legend('$|f(x_k)-f^*|$','${\rm KKT}(x_k)$','${\rm KKT}(y_k)$',...
    '${\rm KKT}(z_k)$','${\rm KKT}(\lambda_k)$'),'Interpreter','latex','FontSize',13);

subplot(1,3,2);
plot(1:length(SsN_itnum),SsN_itnum','hb-','LineWidth',1.5);

xlabel('$k$','Interpreter', 'latex','FontSize',14);
ylabel('$\#{\rm SsN}$','Interpreter', 'latex','FontSize',14);
set(legend('${\rm SsN}(k)$'),'Interpreter','latex','FontSize',14);

subplot(1,3,3);
plot(1:length(PCG_itnum(:,3)),PCG_itnum(:,3)','or-','LineWidth',1.5);hold on;
plot(1:length(PCG_itnum(:,2)),PCG_itnum(:,2)','^b-','LineWidth',1.5);
plot(1:length(PCG_itnum(:,1)),PCG_itnum(:,1)','dk-','LineWidth',1.5);hold off;

xlabel('$k$','Interpreter', 'latex','FontSize',14);
ylabel('$\#{\rm AMG}$','Interpreter', 'latex','FontSize',14);
set(legend('${\rm max\_AMG}(k)$','${\rm ava\_AMG}(k)$','${\rm min\_AMG}(k)$'),...
    'Interpreter','latex','FontSize',14);

figure(2)
plot(1:length(PCG_itnum(:,3)),PCG_itnum(:,3)','or-','LineWidth',1.5);hold on;
plot(1:length(PCG_itnum(:,2)),PCG_itnum(:,2)','^b-','LineWidth',1.5);
plot(1:length(PCG_itnum(:,1)),PCG_itnum(:,1)','dk-','LineWidth',1.5);hold off;

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



