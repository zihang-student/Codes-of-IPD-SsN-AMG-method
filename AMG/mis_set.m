function [isC,isF,As] = mis_set(A,theta)
% Finding the (approximated) maximal independent set by a multilevel type
% algorithm. This is from the ifem toolbox of Long Chen.

%% Parameters
if ~exist('theta','var'), theta = 0.025; end
N = size(A,1);
isF = false(N,1);       % F: fine node
isC = false(N,1);       % C: coarse node
N0 = min(fix(sqrt(N))+1,25);   % number of the coarest nodes

%% Generate strong connectness matrix
% 1. Nonstandard way
% Dinv = spdiags(1./sqrt(diag(A)),0,N,N);
% Am = Dinv*A*Dinv;       % normalize diagonal of A
% [im,jm,sm] = find(Am);
% idx = (-sm > theta);    % delete weakly connect off-diagonal and diagonal
% As = sparse(im(idx),jm(idx),sm(idx),N,N); % matrix for strong connectness
%%
% The diagonal of Am is 1. The negative off-diagonal measures the
% diffusivity. The positive off-diagonal is filtered out.
% 2. standard way
As = strength(A) >= theta;

%% Compute degree of vertex
deg = sum(spones(As)); % number of strongly connected neighbors
deg = full(deg');
if sum(deg>0) < 0.25*sqrt(N)   % too few connected nodes e.g. A is mass matrix
    isC(ceil(rand(N0,1)*N)) = true; % randomly chose N0 nodes
    isF = ~isC;
    return                    % smoother is a good preconditioner
end
idx = (deg>0);deg(idx) = deg(idx) + 0.1*rand(sum(idx),1); % break the
% equal degree case so that every coarsening can chose around half nodes.
% Otherwise the coarsening could be O(N^2). W-cycle is used to make AMG robust.

%% Find an approximate maximal independent set and put to C set
isF(deg == 0) = true;   % isolated nodes are added into F set
isU = true(N,1);        % U: undecided set
while sum(isC) < N/2 && sum(isU) > N0
    % Mark all undecided nodes
    isS = false(N,1);   % S: selected set, changing in the coarsening
    isS(deg>0) = true;  % deg will be set to zero for nodes in C and F
    S = find(isS);
    
    % Find marked nodes with local maximum degree
    [i,j] = find(triu(As(S,S),1));  % i,j and i<j: edges of subgraph S
    idx = deg(S(i)) >= deg(S(j));   % compare degree of connected vertices
    isS(S(j(idx))) = false;         % remove vertices with smaller degree
    isS(S(i(~idx))) = false;        % if degrees are equal, keep the nodes with smaller index
    isC(isS) = true;                % set selected nodes as coarse nodes
    
    % Remove coarse nodes and neighboring nodes from undecided set
    [i,~] = find(As(:,isC)); %#ok<*NASGU>
    isF(i) = true;        % neighbor of C nodes are F nodes
    isU = ~(isF | isC);
    deg(~isU) = 0;        % remove current C and F from the graph
    
    if sum(isU) <= N0     % add small undecided nodes into C nodes
        isC(isU) = true;
        isU = [];         % to exit the while loop;
    end
end
%% Be careful with isolated nodes!!!
iso = sum(As,2)==0;isC(iso) = true;isF(iso) = false;
end