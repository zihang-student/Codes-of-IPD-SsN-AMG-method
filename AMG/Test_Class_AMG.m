

%% A test on classical AMG for solving the linear algebraic system
%                        Ax = b
% where A = eps*I + A0 and A0 is a connected graph Laplacian.

%=======================================================
% clc;
clear;close all;
%% Load graph
load('poly_coarse.mat'); % This provides elem and node
its = 3;[elem,node] = uniref_trian(elem,node,its);
E = [elem(:,[1,2]);elem(:,[2,3]);elem(:,[1,3])];
E = [min(E(:,1),E(:,2)),max(E(:,1),E(:,2))];
E = unique(E,'rows');s = E(:,1);t = E(:,2);

dof = size(node,1);w = -ones(length(s),1);
A = sparse(s,t,w,dof,dof);A = A + A';
A = A + spdiags(abs(sum(A,2)),0,dof,dof);
% G = graph(s,t);IG = incidence(G);
% edge_num = size(G.Edges,1);node_num = size(G.Nodes,1);
% w = rand(edge_num,1);W = spdiags(w,0,edge_num,edge_num);
% w = ones(edge_num,1);W = spdiags(w,0,edge_num,edge_num);
%% Problem setting
epss = 1e-5;A0 = A;
D = spdiags(rand(dof,1),0,dof,dof);
A = epss*D + A0;b = A * A(:,1);
%=======================================================
%% Class_AMG: simple version without display
amg_options = struct('retol',1e-11,'bigph',0,'maxit',5e1,'theta',1/4,...
    'smoth',5,'cycle','v','isnsp',1,'inter',1,'guess',zeros(size(b)));
tic;[~,it,rel_res,rel_resk] = Class_AMG(A,b,amg_options);
disp(['AMG: it = ', num2str(it),'  rel_res = ',num2str(rel_res,'%4.2e'),...
    '  total time  =  ',num2str(toc,'%02.2f'),'s']);
%% Teogrid iteration
% amg_options.fnode = fix(dof/2);
tic;
[~,it,rel_res,rel_resk2] = twogrid(A,b,amg_options);
disp(['Twogrid: it = ', num2str(it),'  rel_res = ',num2str(rel_res,'%4.2e'),...
    '  total time  =  ',num2str(toc,'%02.2f'),'s']);
%% PCG
pcg_options = struct('retol',1e-11,'maxit',1e4,'precd',2,'guess',[]);
tic;[d,it,res,resk] = PCG(A,b,pcg_options);
disp(['PCG: it = ', num2str(it),'  rel_res = ',num2str(res,'%4.2e'),...
    '  time used ',num2str(toc,'%02.2f'),'s']);
%% Class_AMG: old version with more display
amg_options.disp = [1,0,1];amg_options.time = 0;
Class_AMG_detail(A,b,amg_options);
%% Plot
figure;
semilogy(1:length(resk),resk,'bo-');hold on;
semilogy(1:length(rel_resk),rel_resk,'rd-');
semilogy(1:length(rel_resk2),rel_resk2,'kp-');

xlabel('$k$','interpreter','latex','fontsize',15);
ylabel('$\ln \|r_k\|/\|r_0\|$','interpreter','latex','fontsize',15);

set(legend('PCG','AMG','Two-grid'),'Interpreter','latex','FontSize',15);