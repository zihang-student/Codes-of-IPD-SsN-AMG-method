
function y = Ax(x,p,q)
% This code computes the matrix-vector
% multiplication A*x in a matrix-free way, where
%                [ In  ox p']  
%            A = [          ]  
%                [ q' ox Im ]  
% with p in R^m and q in R^n.

m = length(p);n = length(q);
X = reshape(x,m,n);
r = X'*p;l = X*q;
y = [r;l];
end