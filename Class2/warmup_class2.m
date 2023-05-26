%% Warm-start for problem class 2 by using A-ADMM--New version matrix free
function [uk,lk] = warmup_class2(c,r,l,p,q,mu,phi,res,maxit)

if nargin == 7
    res = 1e-1;maxit = inf;
end
if nargin == 8
    maxit = inf;
end
if nargin == 9
    if res == 0 && maxit == inf
        error('res = 0 and maxit = inf');
    end
    if res > 0 && maxit < inf
        warning('require res > 0 and maxit < inf but I chose maxit < inf');
        res = 0;
    end
end
% Preparations
m = length(l);n = length(r);

prox = @(x) max(0,x);b = [r;l;mu];
Htb = [Aty(b(1:m+n),p,q)+b(end)*phi;b(1:n+m)];
wc = [c;zeros(n+m,1)];z0 = zeros(m+n+1,1);

muf = 0; gk = 1;bk = 1;
% Initial guess
uk = spalloc(m*n+n+m,1,0); vk = uk;wk = uk;pik = wk;lk = [z0;uk];
% uk = [c;r;l]; vk = uk;wk = uk;pik = wk;lk = [z0;uk];
% % Residuals
% if maxit ~= inf
%     KKT_lk = - inf*ones(maxit+1,1);
% end
xk = uk(1:m*n);yk = uk(m*n+1:m*n+n);zk = uk(m*n+n+1:end);
% KKT_lk(1) = norm([Ax(xk,p,q)+[yk;zk];phi'*xk]-b);
% KKT_zk = KKT_lk;KKT_zk(1) = norm(zk-max(zk-lk(n+1:n+m),0));
% KKT_yk = KKT_lk;KKT_yk(1) = norm(yk-max(yk-lk(1:n),0));
% KKT_xk = KKT_lk;KKT_xk(1) = norm(xk-max(xk-c-(Aty(lk(1:m+n),p,q)+lk(m+n+1)*phi),0));
% %% 1. Main Loop
% res0 = max(max(max(KKT_xk(1),KKT_yk(1)),KKT_zk(1)),KKT_lk(1));tic;
disp(['===========Warm-start===== ',datestr(now)]);
% aa = 'KKT(xk) = %4.2e, KKT(yk) = %4.2e, KKT(zk) = %4.2e, KKT(lk) = %4.2e\n';
% dd1 = ['A-ADMM: it = %4d, ',aa];dd2 = ['        it = %4d, ',aa];
% tic;
k = 1;
while k <= maxit
    %     resk = max(max(max(KKT_xk(k),KKT_yk(k)),KKT_zk(k)),KKT_lk(k));
    %     if k == 1 || k == maxit
    %         fprintf(dd1,k-1,KKT_xk(k),KKT_yk(k),KKT_zk(k),KKT_lk(k));
    %     else
    %         if res0/resk >= 2
    %             res0 = resk;
    %             fprintf(dd2,k-1,KKT_xk(k),KKT_yk(k),KKT_zk(k),KKT_lk(k));
    %         end
    %     end
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
    % 1.3 Traansfer
    gk = gk1;bk = bk1;uk = uk1;vk = vk1;
    wk = wk1;pik = pik1;lk = lk1;
    xk = uk(1:m*n);yk = uk(m*n+1:m*n+n);zk = uk(m*n+n+1:end);
    
    if mod(k,50) == 0
        gg = 'A-ADMM it = %3d\n';fprintf(gg,k);
    end
    
    %     KKT_lk(k+1) = norm([Ax(xk,p,q)+[yk;zk];phi'*xk]-b);
    %     KKT_zk(k+1) = norm(zk-max(zk-lk(n+1:n+m),0));
    %     KKT_yk(k+1) = norm(yk-max(yk-lk(1:n),0));
    %     KKT_xk(k+1) = norm(xk-max(xk-c-(Aty(lk(1:m+n),p,q)+lk(m+n+1)*phi),0));
    %
    %     rr1 = [KKT_xk(k+1)/(1+KKT_xk(1)),KKT_yk(k+1)/(1+KKT_yk(1))];
    %     rr = [rr1 KKT_zk(k+1)/(1+KKT_zk(1)),KKT_lk(k+1)/(1+KKT_lk(1))];
    %     if max(rr) <= res
    %         gg = 'A-ADMM stopped at it = %3d, rel_KKT_res = %4.2e, time used %4.2es\n';
    %         fprintf(gg,k,max(rr),toc);break;
    %     end
    %     if k == maxit
    %         gg = 'A-ADMM stopped at MAXit = %3d, rel_KKT_res = %4.2e, time used %4.2es\n';
    %         fprintf(gg,k,max(rr),toc);break;
    %     end
    k = k + 1;
end
lk = lk(1:m+n+1);
end