
function [e] = MG_Vcycle(r,isnsp,k)
%% Recursive MG-Vcycle for solving Ae = r with zero initial guess

if nargin == 1;isnsp = 0;k = 1;end
if nargin == 2;k = 1;end

global Ack Prok J smoth_it Rk;
R = Rk{k};A = Ack{k};
N = size(A,1);Rt = R';

if k < J
    % Presmoothing
    if isnsp % add kernel space
        xi = ones(N,1);xx = xi'*A*xi;Axi = A*xi;
        e = zeros(size(r));
        for i = 1:smoth_it
            g = r-A*e;xig = xi'*g;
            g = xi*(xig/xx) + R*(g-Axi*(xig/xx));
            e = e + g;
        end
    else
        e = zeros(size(r));
        for i = 1:smoth_it; e = e + R*(r-A*e);end
    end
    % Restriction
    rr = r-A*e;rrc = Prok{k+1}'*rr;
    % Correction
    eec = MG_Vcycle(rrc,isnsp,k+1);
    % Prolongation
    ee = Prok{k+1}*eec;e = e + ee;
    % Postmoothing
    if isnsp  % add kernel space
        for i = 1:smoth_it
            g = r-A*e;xig = xi'*g;
            g = xi*(xig/xx) + Rt*(g-Axi*(xig/xx));
            e = e + g;
        end
    else
        for i = 1:smoth_it;e = e + Rt*(r-A*e);end
    end
else % k is the coasest level
    [e,~] = PCG(A,r);
    %     e = A \ r;
end
end