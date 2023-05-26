function [indC,indF,SubG] = cf_split(S)
%% This returns the C/F splitting information based on
%   the given logical strength connection matrix S.

%% Build the subgraph w.r.t. the strength connection matrix S
SubG = graph(S);
N = size(S,1);indF = false(N,1);
indC = false(N,1);indU = true(N,1);

for k = 1 : N
    kk = neighbors(SubG,k);
    if indU(k) == true % node k has not been visited
        indC(k) = true;indU(k) = false;
        indF(kk) = true;indU(kk) = false;
    end
end