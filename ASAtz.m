
function y = ASAtz(z,s,p,q)
% This code computes the matrix-vector multiplication 
% A*diag(s)*A'*z in a matrix-free way, where
%                [ In  ox p']  
%            A = [          ]  
%                [ q' ox Im ]  
% with p in R^m, q in R^n and s in {0,1}^{m*n}.

% Let vec(Y) = s, U = diag(p)*Y and Q = Y*diag(q), then 
%                 [ diag(U'*p)  diag(q)*U']  
%  A*diag(s)*A' = [                       ]  
%                 [ diag(p)*Q   diag(Q*p) ] 

m = length(p);n = length(q);
Y = reshape(s,m,n);
U = diag(p)*Y;Q = Y*diag(q);
z1 = z(1:n);z2 = z(n+1:n+m);

y1 = (U'*p).*z1 + q.*(U'*z2);
y2 = p.*(Q*z1) + (Q*p).*z2;
y = [y1;y2];
end
