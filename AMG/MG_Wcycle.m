
function [e] = MG_Wcycle(r,isnsp,k,e)
%% Recursive MG-Wcycle for solving Ae = r

if nargin == 1;isnsp = 0;k = 1;e = zeros(size(r));end
if nargin == 2;k = 1;e = zeros(size(r));end
if nargin == 3;e = zeros(size(r));end

global Ack Prok J smoth_it Rk;
R = Rk{k};A = Ack{k};
N = size(A,1);Rt = R';

if k < J
    % Presmoothing
    if isnsp
        xi = ones(N,1);xx = xi'*A*xi;Axi = A*xi;
        for i = 1:smoth_it
            g = r-A*e;xig = xi'*g;
            g = xi*(xig/xx) + R*(g-Axi*(xig/xx));
            e = e + g;
        end
    else
        for i = 1:smoth_it; e = e + R*(r-A*e);end
    end
    % Restriction
    rr = r-A*e;rrc = Prok{k+1}'*rr;
    % Correction
    eec = MG_Wcycle(rrc,isnsp,k+1);
    % Correction again
    eec = MG_Wcycle(rrc,isnsp,k+1,eec);
    % Prolongation
    ee = Prok{k+1}*eec;e = e + ee;
    % Postmoothing
    if isnsp
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
    % e = A \ r;
end
end