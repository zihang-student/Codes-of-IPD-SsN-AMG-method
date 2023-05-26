function y = invHHt(v,p,q,sg,phi)
% This code computes y = (sg*I + H*H') \ x where H = (G,IY,IZ), 
% G = [A;phi'] and
%                 [ In  ox p']
%             A = [          ]
%                 [ q' ox Im ]
m = length(p);n = length(q);
t = sg + norm(phi)^2;l = Ax(phi,p,q);
Vl = invAAt(l,p,q,sg+1);s = t - l'*Vl;
%%
v1 = v(1:n+m);v2 = v(end);
Vv1 = invAAt(v1,p,q,sg+1);
%%
y1 = s*Vv1 + l'*Vv1 * Vl - v2*Vl;
y2 = v2 - l'*Vv1;
%%
y = [y1;y2]/s;
end
