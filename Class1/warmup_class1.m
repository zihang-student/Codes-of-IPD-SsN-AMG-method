%% Warm-start for problem class 1 by using A-ADMM--matrix free
function [xk,lk] = warmup_class1(c,r,l,p,q,gama,res,maxit)
if nargin == 6
    res = 1e-1;maxit = inf;
end
if nargin == 7
    maxit = inf;
end
if nargin == 8
    if res == 0 && maxit == inf
        error('res = 0 and maxit = inf');
    end
    if res > 0 && maxit < inf
        warning('require res > 0 and maxit < inf but I chose maxit < inf');
        res = 0;
    end
end
if maxit == inf
    maxit = 500;
end
% Preparations
m = length(l);n = length(r);

prox = @(x) min(max(0,x),gama);
b = [r;l];Atb = Aty(b,p,q);z0 = zeros(m+n,1);

muf = 0; gk = 1;bk = 1;
% Initial guess
xk = spalloc(m*n,1,0);vk = xk;wk = xk;pik = wk;lk = [z0;xk];
% Residuals
% if maxit ~= inf
%     KKT_lk = - inf*ones(maxit+1,1);
% end
% KKT_lk(1) = norm(Ax(xk,p,q)-b);KKT_xk = KKT_lk;
% KKT_xk(1) = norm(xk-prox(xk-c-Aty(lk(1:m+n),p,q)));
%% 1. Main Loop
% res0 = max(KKT_xk(1),KKT_lk(1));tic;
disp(['===========Warm-start===== ',datestr(now)]);
% aa = 'A-ADMM: it = %6d, KKT(xk) = %4.2e, KKT(lk) = %4.2e\n';
% bb = '        it = %6d, KKT(xk) = %4.2e, KKT(lk) = %4.2e\n';
% tic;
k = 1;
while k <= maxit
    %     resk = max(KKT_xk(k),KKT_lk(k));
    %     if k == 1 || k == maxit
    %         if k == 1
    %             fprintf(aa,k-1,KKT_xk(k),KKT_lk(k));
    %         else
    %             fprintf(bb,k-1,KKT_xk(k),KKT_lk(k));
    %         end
    %     else
    %         if res0/resk >= 2
    %             res0 = resk;fprintf(bb,k-1,KKT_xk(k),KKT_lk(k));
    %         end
    %     end
    % 1.1 Iterate
    ak = bk;bk1 = bk/(1+ak);
    gk1 = (gk+muf*ak)/(1+ak);
    etafk = (1+ak)*gk + muf*ak;
    sgk = 1/bk1;etagk = (1+ak)*bk;
    
    wwk = (ak*pik + wk)/(1+ak);
    wxk = (ak*gk*vk + (gk+muf*ak)*xk)/etafk;
    
    hlk = lk - 1/bk*[Ax(xk,p,q)-b;xk-wk] + ak/bk*[z0;-(pik-wk)];
    cAw = - Atb - wk;cAlk = Aty(hlk(1:m+n),p,q)+hlk(m+n+1:end);
    dd = etafk*wxk - ak^2*(c+cAlk+sgk*cAw);
    % 1.2 Update
    tt = sgk*ak^2;sg = 1+etafk/tt;
    xk1 = (dd-Aty(invAAt(Ax(dd,p,q),p,q,sg),p,q))/(etafk+tt);
    vk1 = xk1 + (xk1 - xk)/ak;
    blk = lk + ak/bk*[Ax(vk1,p,q)-b;vk1-pik];
    wk1 = prox(wwk-ak^2/etagk*(-blk(m+n+1:end)));
    pik1 = wk1 + (wk1-wk)/ak;
    lk1 = lk + ak/bk*[Ax(vk1,p,q)-b;vk1-pik1];
    % 1.3 Transfer
    gk = gk1;bk = bk1;xk = xk1;vk = vk1;wk = wk1;pik = pik1;lk = lk1;
    
    if mod(k,50) == 0
        gg = 'A-ADMM it = %3d\n';fprintf(gg,k);
    end
    %     KKT_lk(k+1) = norm(Ax(xk,p,q)-b);
    %     KKT_xk(k+1) = norm(xk-prox(xk-c-Aty(lk(1:m+n),p,q)));
    %
    %     rr = [KKT_xk(k+1)/(1+KKT_xk(1)),KKT_lk(k+1)/(1+KKT_lk(1))];
    %     if max(rr) <= res
    %         gg = 'A-ADMM stopped at it = %3d, KKT_res = %4.2e, time used %4.2es\n';
    %         fprintf(gg,k,max(rr),toc);break;
    %     end
    %     if k == maxit
    %         gg = 'A-ADMM stopped at MAXit = %3d, KKT_res = %4.2e, time used %4.2es\n';
    %         fprintf(gg,k,max(rr),toc);break;
    %     end
    k = k + 1;
end
lk = lk(1:m+n);
end