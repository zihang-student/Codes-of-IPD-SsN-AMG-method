function [zeta,itpcg,respcg,info] = PCG4POT(prob_data,pcg_options)
% This is only for partial optimal transport (POT).

% Original system:
%                   He * zeta = z
% with He = eps*I + sg*H and
%              [T + A*S*A'    A*S*phi  ]
%           H =[                       ]
%              [phi'*S*A'    phi'*S*phi]
% where S = diag(s) and T = diag(t).

% Let z = (z1,z2), phi_e = eps + sg*phi'*S*phi and w = z1 - sg/phi_e * z2*v.
% It is equivalent to
%     zeta1 = (Ae - sg^2/phi_e * v*v') \ w,
%     zeta2 = (z2 - sg*v' * zeta1)/phi_e,
% where v = A*S*phi and Ae = eps*I + sg*(T + A*S*A').

% By Sherman-Woodbury formula, we have
%   inv(Ae - sg^2/phi_e * v*v') = inv(Ae) + tt * inv(Ae) * v*v' * inv(Ae),
% where tt = sg^2/(phi_e - sg^2 * v'*inv(Ae)*v).
% Then we have to solve two linear systems
%          Ae * ww = w and Ae * vv = v.
% This gives zeta1 = ww + tt * vv*v'*ww with
%        tt = sg^2/(phi_e - sg^2 * v'*vv).

p = prob_data.p; q = prob_data.q;
bk1 = prob_data.bk1;tk = prob_data.tk;
phi = prob_data.phi;z = prob_data.z;s = prob_data.s;

z1 = z(1:end-1);z2 = z(end); epss = bk1; sg = 1/tk;

phi_e = epss + sg * phi'*(s.*phi);
v = Ax(s.*phi,p,q);w = z1 - sg/phi_e * z2*v;

prob_data.z = v;[vv,itpcg1,respcg1,info1] = aug_PCG(prob_data,pcg_options);
prob_data.z = w;[ww,itpcg2,respcg2,info2] = aug_PCG(prob_data,pcg_options);

tt = sg^2/(phi_e - sg^2*v'*vv);
zeta1 = ww + tt*vv*v'*ww;zeta2 = (z2-sg*v'*zeta1)/phi_e;zeta = [zeta1;zeta2];
itpcg = max(itpcg1,itpcg2); respcg = max(respcg1,respcg2);info = max(info1,info2);
end