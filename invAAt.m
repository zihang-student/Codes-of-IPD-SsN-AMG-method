function y = invAAt(x,p,q,sg1,sg2)
% This code computes y = (diag(sg1,sg2) + A*A') \ x where
%              [ In  ox p']
%          A = [          ]
%              [ q' ox Im ]
%!!!!!!! This is not robust w.r.t. sg.
if nargin == 3
    sg1 = 1;sg2 = 1;
end
if nargin == 4
    sg2 = sg1;
end
m = length(p);np = norm(p)^2;n = length(q);nq = norm(q)^2;
%%
vn = x(1:n);vm = x(n+1:n+m);
%%
yn = vn/(sg1+np) + (np/(sg1+np)*(q'*vn)-p'*vm)*q/(sg1*sg2+sg1*nq+sg2*np);
ym = vm/(sg2+nq) + (nq/(sg2+nq)*(p'*vm)-q'*vn)*p/(sg1*sg2+sg1*nq+sg2*np);
%%
y = [yn;ym];
end
