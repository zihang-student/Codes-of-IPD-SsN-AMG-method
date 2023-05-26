function [zeta,itpcg,respcg,info] = aug_PCG(prob_data,pcg_options)
% Original system:
%                   He * zeta = z
% with He = bk1*I + 1/tk*H, H = T + H0, and
%         S = diag(s), T = diag(t), H0 = A*S*A'.
% Let zeta = Q0 * u, we have an equivalent system:
%                     Ae * u = f
% with Ae = Q0*He*Q0 = bk1*Q0^2 + 1/tk*(K+A0) and f = Q0*z, where
%         Q0 = diag(q,-p), K = Q0*T*Q0, A0 = Q0*H0*Q0.
%% He * zeta = z
bk1 = prob_data.bk1;tk = prob_data.tk;
q = prob_data.q;p = prob_data.p;
H0 = prob_data.H0;z = prob_data.z;
T = prob_data.T;M = size(H0,1);
%% Ae * u = f
qp = [q;-p];
if any(qp==0)
    error('p or q contains 0 !!!!!');
else
    Q0 = spdiags(qp,0,M,M);
end
A0 = Q0*H0*Q0; Q = Q0*Q0; K = Q0*T*Q0; f = Q0*z;
%% Null space of A0
[blocks,sizes] = components(A0); % This is from the Lean-AMG toolbox
num_comp = length(sizes);Y = sparse(1:M,blocks,ones(M,1),M,num_comp);
%% Solve the augmented system via PCG
Ae = bk1*Q + 1/tk*(K+A0);QK = bk1*Q + 1/tk*K; % A0*Y = 0
augAe = [(Y')*QK*Y Y'*QK;QK*Y Ae];augf = [Y'*f;f];
pcg_options.guess = [zeros(size(Y,2),1);zeros(size(f))];

% warning('bi-SSOR has not been implemented for augmented system!!');
pcg_options.precd = 2;

[U,itpcg,respcg] = PCG(augAe,augf,pcg_options);
u = Y*U(1:size(Y,2)) + U(size(Y,2)+1:end);
%% zeta
zeta = Q0*u;info = [num_comp,1];
end