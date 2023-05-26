
function z = Aty(y,p,q)
% This code computes the matrix-vector
% multiplication A'*y in a matrix-free way, where
%                [ In  ox p']  
%            A = [          ]  
%                [ q' ox Im ]  
% with p in R^m and q in R^n.

m  = length(p);n = length(q);
y1 = y(1:n);y2 = y(n+1:n+m);
z1 = p*y1';z2 = y2*q';
z = z1(:) + z2(:);
end