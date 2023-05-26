function [S] = strength(A,which)
%% This builds the strength value matrix S based on the given matrix A.

if nargin == 1;which = 2;end

N = size(A,1);
D = spdiags(diag(A),0,N,N);
A0 = D-A;[ia,ja,s0] = find(A0);
max_row = max(A0,[],2);
max_row(max_row<=0) = inf;

switch which
    case 1  %% This is not symmetric!!!!!
        sa = s0./max_row(ia);
    case 2  %% Symmetrized version of case 1
        sa = s0./min(max_row(ia),max_row(ja)); 
end
S = sparse(ia,ja,sa,N,N); % not logical but strength value !!!!!!
end