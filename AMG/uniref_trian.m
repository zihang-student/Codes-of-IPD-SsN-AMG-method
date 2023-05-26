function [elem,node] = uniref_trian(elem,node,its)
% Refine a given triangular mesh by uniform refinement with its times.

if nargin < 3
    its = 1; % Only one refinement
end

while its
    node_num = size(node,1);elem_num = size(elem,1);
    edge     = sort([elem(:,1:2);elem(:,2:3);elem(:,[1,3])],2);
    
    [edge,~,elem2edge] = unique(edge,'rows') ;
    elem2edge = reshape(elem2edge,elem_num,3);
    
    % % Elements t = [P1,P2,P3], edges te = [P1-P2,P2-P3 P3-P1]
    
    %                 P3
    %                /  \
    %               /    \
    %              /  3   \
    %           e3/________\e2
    %            / \      / \
    %           /   \  4 /   \
    %          / 1   \  / 2   \
    %         /_______\/_______\
    %        P1       e1        P2
    
    % New element
    new_elem = zeros(4*elem_num,3);
    % P1  Pe1 Pe3
    new_elem(1:4:4*elem_num-3,:) = ...
        [elem(:,1),node_num+[elem2edge(:,1),elem2edge(:,3)]];
    % P2  Pe2 Pe1
    new_elem(2:4:4*elem_num-2,:) = ...
        [elem(:,2),node_num+[elem2edge(:,2),elem2edge(:,1)]];
    % P3  Pe3 Pe2
    new_elem(3:4:4*elem_num-1,:) = ...
        [elem(:,3),node_num+[elem2edge(:,3),elem2edge(:,2)]];
    % Pe1 Pe2 Pe3
    new_elem(4:4:4*elem_num,  :) = ...
        [elem2edge(:,1),elem2edge(:,2),elem2edge(:,3)] + node_num;
    elem = new_elem;
    % New node
    new_node = [node;(node(edge(:,1),:) + node(edge(:,2),:))/2];
    node = new_node;its = its - 1;
end
end