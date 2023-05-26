function [d,it,res,resk] = PCG(H,e,pcg_options)
% This is the PCG method for solving the linear SPD system Hd = e.
% It can be found in Appendix B3 of the paper
% ﻿J. Shewchuk. An introduction to the conjugate gradient method without the
%  agonizing, edition 5/4. Technical report, Carnegie Mellon University,
%  Pittsburgh, PA, USA, 1994.
%%====================================
% pcg_tol: pcg_tolerance of relative residual
% pcg_maxit: maximal iteration number
% d0: initial guess
% ii: choice of preconditioner
%  ii = 1 P = I: No preconditioner
%  ii = 2 P = D: Jacobi preditioner
%  ii = 3 P = SSOR: Symmetric SOR
%  ii = 4 P = iC: incomplete Choleskey
%  ii = 5 P = bi-SSOR: SSOR only for bigraph Laplacian
%%====================================
if nargin == 2
    pcg_options.guess = zeros(size(e));
    pcg_options.retol = 1e-11;
    pcg_options.maxit = 1e4;
    pcg_options.precd = 2;
end
if isempty(pcg_options.guess); pcg_options.guess = zeros(size(e));end
if isempty(pcg_options.retol); pcg_options.retol = 1e-11;end
if isempty(pcg_options.maxit); pcg_options.maxit = 1e4;end
if isempty(pcg_options.precd); pcg_options.precd = 2;end

ii = pcg_options.precd;
d0 = pcg_options.guess;
pcg_tol = pcg_options.retol;
pcg_maxit = pcg_options.maxit;

switch ii
    case 1 % No preconditioner
        P = 1;
    case 2 % Jacobi preditioner
        P = reshape(diag(H),size(d0));
    case 3 % Symmetric SOR
        P{2} = tril(H,-1);P{3} = triu(H,1);
        P{1} = spdiags(diag(H),0,size(H,1),size(H,2));
        % P = 1/(w*(2-w)) * (D+w*L)*inv(D)*(D+w*U) 1<= w < 2
        % P^{-1} = (w*(2-w)) * (D+w*U)^{-1}*D*(D+w*L)^{-1}
    case 4 % incomplete Choleskey
        if issparse(H)
            P = ichol(H);
            %             P = ichol(H,struct('type','ict','droptol',1e-02,'michol','on'));
        else
            error('iC requires H is sparse!');
        end
    case 5 % Symmetric SOR for bigraph Laplacian
        %     [V    U]
        % H = [      ]
        %     [U'   T]
        if isfield(pcg_options,'nf')
            w = 1.5;Nf = pcg_options.nf;
            V = H(1:Nf,1:Nf);U = H(1:Nf,Nf+1:end);
            T = H(Nf+1:end,Nf+1:end);Nc = size(H,1)-Nf;
            invV = spdiags(1./diag(V),0,Nf,Nf);
            invT = spdiags(1./diag(T),0,Nc,Nc);
            P = w*(2-w)*[invV + w^2*invV*U*invT*U'*invV,...
                -w*invV*U*invT;-w*invT*U'*invV, invT];
        else
            error('SSOR for bigraph requires pcg_options.nf!!!');
        end
end

it = 0; r = e - H*d0;
p = pre_cond_M(r,ii,P);
delta_new = r'*p;d = d0;

delta_0 = delta_new;
% delta_0 = 1;
resk = zeros(pcg_maxit,1);
% Be careful with the stagnation
while it < pcg_maxit && delta_new > pcg_tol^2*delta_0
    delta_old = delta_new;q = H*p;
    alpha = delta_old/(q'*p);d = d+alpha*p;
    r = r - alpha*q;
    w = pre_cond_M(r,ii,P);
    delta_new = r'*w;
    beta = delta_new/delta_old;
    p = w + beta*p;it = it + 1;
    
    resk(it) = sqrt(abs(delta_new/delta_0));
end
res = sqrt(abs(delta_new/delta_0));
end

function p = pre_cond_M(r,ii,P)
switch ii
    case 1
        p = r;
    case 2
        p = r./P;
    case 3
        w = 1.5;
        p1 = (P{1}+w*P{2})\r;p2 = P{1}*p1;
        p = w*(2-w) * (P{1}+w*P{3}) \ p2;
    case 4
        p = P\r;p = P'\p;
    case 5
        p = P*r;
end
end