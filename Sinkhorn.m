%% This code tests the Sinkhorn.
% X = exp(u/eps) * exp(v'/eps) .* exp(-1-C/eps)

% clc;clear;close all;
C = reshape(c,m,n);
% ic = 2; 
% n = 1e3;
tic;
%% 1. Chose the cost matrix
% % 1. random cost
% if ic == 1
%     C = rand(n,n);c = abs(C(:));
% end
% % 2. {0,1}-cost
% if ic == 2
%     C = ones(n) - eye(n); c = abs(C(:)); % half l1-norm = 0.5*|b1-b2|_1
% end
% % 3. lp-cost
% if ic == 3
%     dim = 10;p = 1;
%     x = rand(dim,1,n);ii = (1:n)'*ones(1,n);
%     y = rand(dim,1,n);jj =  ones(n,1)*(1:n);
%     ss = abs(x(:,:,ii)-y(:,:,jj));
%     ss = sum(ss.^p);ss=ss(:);
%     C = sparse(ii,jj,ss,n,n);
%     c = abs(C(:));
% end
% %% 2. Generate the source vector a and the target vector b
% d = rand(2*n,1);
% if sum(d(1:n))>sum(d(n+1:end))
%     d = [d(1:n);d(n+1:end)+(sum(d(1:n))-sum(d(n+1:end)))/n];
% else
%     d = [d(1:n)-(sum(d(1:n))-sum(d(n+1:end)))/n;d(n+1:end)];
% end
% l = d(1:n);r = d(n+1:end);
%% 3. Sinkhorn 
eps = 1e-3;Keps = exp(-1-C/eps);

uk = -rand(n,1);vk = -rand(n,1);
dk = uk'*l+vk'*r -eps*exp(uk'/eps)*Keps*exp(vk/eps);

K = 1e3;k = 0;
if ic ~= 2
    Kk = 2*K;
else
    Kk = K;
end

while k < Kk
    uk1 = eps * log(l./(Keps *exp(vk /eps)));
    vk1 = eps * log(r./(Keps'*exp(uk1/eps)));    
    
    k = k+1;uk = uk1;vk = vk1;
    
    dkk = uk'*l+vk'*r -eps*exp(uk'/eps)*Keps*exp(vk/eps);
    dk0 = dk;dk = [dk0,dkk];
end
X = exp((uk+vk'-C)/eps-1); X(X<1e-3) = 0;
toc;