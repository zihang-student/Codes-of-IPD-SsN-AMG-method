
function H = ASAt(s,p,q)
% This returns the (sparse) matrix A*diag(s)*A', where
%                [ In  ox p']
%            A = [          ]
%                [ q' ox Im ]
% with p in R^m, q in R^n and s in {0,1}^{m*n}.

% Let vec(Y) = s, U = diag(p)*Y and Q = Y*diag(q), then
%                 [ diag(U'*p)  diag(q)*U']
%  A*diag(s)*A' = [                       ]
%                 [ diag(p)*Q   diag(Q*q) ]

m = length(p);n = length(q);
Y = sparse(reshape(s,m,n));
P = spdiags(p,0,m,m);
R = spdiags(q,0,n,n);
U = P*Y;Q = Y*R;
H = [spdiags(U'*p,0,n,n) R*U';P*Q spdiags(Q*q,0,m,m)];
end
