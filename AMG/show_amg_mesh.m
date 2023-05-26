function  show_amg_mesh(G,S,indC,F2C,detail)

if nargin == 4
    detail = false;
end
figure;J = length(G);
%% Visualization
if detail
    %% Level = 1
    subplot(2,4,1);
    Gk = G{1}; Gk =  Gk - spdiags(diag(Gk),0,size(Gk,1),size(Gk,1));
    Gk = tril(Gk);Gk = graph(Gk + Gk');
    indCk = indC{1};Sk = graph(S{1});
    % connected submesh by S in original mesh by G
    [~,ind1]= intersect(Gk.Edges.EndNodes,Sk.Edges.EndNodes,'rows');
    N_edge = size(Gk.Edges.EndNodes,1);
    EdgeCData = ones(N_edge,1)*[0 0 1]; % blue for all edge
    EdgeCData(ind1,:) = ones(size(ind1,1),1)*[0 0 0]; % for strong connected edge
    wt = ones(N_edge,1); wt(ind1) = 3;
    
    % Node color
    N_node = size(Gk.Nodes,1);
    C_node = find(indCk);num_c = length(C_node);
    NodeCData = ones(N_node,1)*[0 0 1];
    NodeCData(C_node,:) = ones(num_c,1)*[1 0 0]; % red for Cnode
    % Node mark size
    NodeMData = 2.5*ones(N_node,1);
    NodeMData(C_node,:) = 3.5*ones(num_c,1);
    
    
    plot(Gk,'MarkerSize',NodeMData,'NodeColor',NodeCData,...
        'EdgeColor',EdgeCData,'Layout',...
        'force','Iterations',20,'LineWidth',wt);axis off;
    title(['Level = ','$1 ',':~',' G$'],'Interpreter','latex','fontsize',14);
    % connected submesh by S
    subplot(2,4,2);F2Ck = F2C{1};
    F2C_edge = [(1:length(F2Ck))',F2Ck];
    F2C_edge(F2C_edge(:,2) == 0,:) = [];
    F2C_edge = [min(F2C_edge,[],2),max(F2C_edge,[],2)];
    [~,ind2]= intersect(Sk.Edges.EndNodes,F2C_edge,'rows');
    
    N_edge = size(Sk.Edges.EndNodes,1);
    EdgeCData = ones(N_edge,1)*[0 0 1]; % blue for all edge
    EdgeCData(ind2,:) = ones(size(ind2,1),1)*[0 0 0]; % for strong connected edge
    wt = ones(N_edge,1); wt(ind2) = 3;
    
    plot(Sk,'MarkerSize',NodeMData,'NodeColor',NodeCData,...
        'EdgeColor',EdgeCData,'Layout','force','Iterations',20,'LineWidth',wt);axis off;
    title(['Level = ','$1 ',':~',' S$'],'Interpreter','latex','fontsize',14);
    %% Level = 2
    subplot(2,4,3);
    Gk = G{2}; Gk =  Gk - spdiags(diag(Gk),0,size(Gk,1),size(Gk,1));
    Gk = tril(Gk);Gk = graph(Gk + Gk');
    
    indCk = indC{2};Sk = graph(S{2});
    % connected submesh by S in original mesh by G
    [~,ind1]= intersect(Gk.Edges.EndNodes,Sk.Edges.EndNodes,'rows');
    N_edge = size(Gk.Edges.EndNodes,1);
    EdgeCData = ones(N_edge,1)*[0 0 1]; % blue for all edge
    EdgeCData(ind1,:) = ones(size(ind1,1),1)*[0 0 0]; % for strong connected edge
    wt = ones(N_edge,1); wt(ind1) = 3;
    
    % Node color
    N_node = size(Gk.Nodes,1);
    C_node = find(indCk);num_c = length(C_node);
    NodeCData = ones(N_node,1)*[0 0 1];
    NodeCData(C_node,:) = ones(num_c,1)*[1 0 0]; % red for Cnode
    % Node mark size
    NodeMData = 2.5*ones(N_node,1);
    NodeMData(C_node,:) = 3.5*ones(num_c,1);
    
    
    plot(Gk,'MarkerSize',NodeMData,'NodeColor',NodeCData,...
        'EdgeColor',EdgeCData,'Layout',...
        'force','Iterations',20,'LineWidth',wt);axis off;
    title(['Level = ','$2 ',':~',' G$'],'Interpreter','latex','fontsize',14);
    % connected submesh by S
    subplot(2,4,4);F2Ck = F2C{2};
    F2C_edge = [(1:length(F2Ck))',F2Ck];
    F2C_edge(F2C_edge(:,2) == 0,:) = [];
    F2C_edge = [min(F2C_edge,[],2),max(F2C_edge,[],2)];
    [~,ind2]= intersect(Sk.Edges.EndNodes,F2C_edge,'rows');
    
    N_edge = size(Sk.Edges.EndNodes,1);
    EdgeCData = ones(N_edge,1)*[0 0 1]; % blue for all edge
    EdgeCData(ind2,:) = ones(size(ind2,1),1)*[0 0 0]; % for strong connected edge
    wt = ones(N_edge,1); wt(ind2) = 3;
    
    plot(Sk,'MarkerSize',NodeMData,'NodeColor',NodeCData,...
        'EdgeColor',EdgeCData,'Layout','force','Iterations',20,'LineWidth',wt);axis off;
    title(['Level = ','$2 ',':~',' S$'],'Interpreter','latex','fontsize',14);
    if J == 3
        %% Level = J
        subplot(2,4,5);Gk = G{J};
        Gk =  Gk - spdiags(diag(Gk),0,size(Gk,1),size(Gk,1));
        Gk = tril(Gk);Gk = graph(Gk + Gk');
        
        plot(Gk,'Layout','force','Iterations',50);axis off;
        title(['Level = ','$J',' = ',num2str(J),':~',' G$'],'Interpreter','latex','fontsize',14);
    end
    if J > 3
        %% Level = J-1
        subplot(2,4,5);
        Gk = G{J-1};
        Gk =  Gk - spdiags(diag(Gk),0,size(Gk,1),size(Gk,1));
        Gk = tril(Gk);Gk = graph(Gk + Gk');
        
        indCk = indC{J-1};Sk = graph(S{J-1});
        % connected submesh by S in original mesh by G
        [~,ind1]= intersect(Gk.Edges.EndNodes,Sk.Edges.EndNodes,'rows');
        N_edge = size(Gk.Edges.EndNodes,1);
        EdgeCData = ones(N_edge,1)*[0 0 1]; % blue for all edge
        EdgeCData(ind1,:) = ones(size(ind1,1),1)*[0 0 0]; % for strong connected edge
        wt = ones(N_edge,1); wt(ind1) = 3;
        
        % Node color
        N_node = size(Gk.Nodes,1);
        C_node = find(indCk);num_c = length(C_node);
        NodeCData = ones(N_node,1)*[0 0 1];
        NodeCData(C_node,:) = ones(num_c,1)*[1 0 0]; % red for Cnode
        % Node mark size
        NodeMData = 2.5*ones(N_node,1);
        NodeMData(C_node,:) = 3.5*ones(num_c,1);
        
        
        plot(Gk,'MarkerSize',NodeMData,'NodeColor',NodeCData,...
            'EdgeColor',EdgeCData,'Layout',...
            'force','Iterations',20,'LineWidth',wt);axis off;
        title(['Level = ','$J-1',' = ',num2str(3),':~',' G$'],'Interpreter','latex','fontsize',14);
        % connected submesh by S
        subplot(2,4,6);F2Ck = F2C{J-1};
        F2C_edge = [(1:length(F2Ck))',F2Ck];
        F2C_edge(F2C_edge(:,2) == 0,:) = [];
        F2C_edge = [min(F2C_edge,[],2),max(F2C_edge,[],2)];
        [~,ind2]= intersect(Sk.Edges.EndNodes,F2C_edge,'rows');
        
        N_edge = size(Sk.Edges.EndNodes,1);
        EdgeCData = ones(N_edge,1)*[0 0 1]; % blue for all edge
        EdgeCData(ind2,:) = ones(size(ind2,1),1)*[0 0 0]; % for strong connected edge
        wt = ones(N_edge,1); wt(ind2) = 3;
        
        plot(Sk,'MarkerSize',NodeMData,'NodeColor',NodeCData,...
            'EdgeColor',EdgeCData,'Layout','force','Iterations',20,'LineWidth',wt);axis off;
        title(['Level = ','$J-1',' = ',num2str(J-1),':~',' S$'],'Interpreter','latex','fontsize',14);
        %% Level = J
        subplot(2,4,7);Gk = G{J};
        Gk =  Gk - spdiags(diag(Gk),0,size(Gk,1),size(Gk,1));
        Gk = tril(Gk);Gk = graph(Gk + Gk');
        plot(Gk,'Layout','force','Iterations',50);axis off;
        title(['Level = ','$J',' = ',num2str(J),':~',' G$'],'Interpreter','latex','fontsize',14);
    end
else
    subplot(2,2,1);Gk = G{1};
    Gk =  Gk - spdiags(diag(Gk),0,size(Gk,1),size(Gk,1));
    Gk = tril(Gk);Gk = graph(Gk + Gk');
    plot(Gk,'Layout','force','Iterations',50);axis off;
    title(['Level = ','$1$'],'Interpreter','latex','fontsize',14);
    
    subplot(2,2,2);Gk = G{2};
    Gk =  Gk - spdiags(diag(Gk),0,size(Gk,1),size(Gk,1));
    Gk = tril(Gk);Gk = graph(Gk + Gk');
    plot(Gk,'Layout','force','Iterations',50);axis off;
    title(['Level = ','$2$'],'Interpreter','latex','fontsize',14);
    if J == 3
        subplot(2,2,3);Gk = G{J};
        Gk =  Gk - spdiags(diag(Gk),0,size(Gk,1),size(Gk,1));
        Gk = tril(Gk);Gk = graph(Gk + Gk');
        plot(Gk,'Layout','force','Iterations',50);axis off;
        title(['Level = ','$J$',' = ',num2str(3)],'Interpreter','latex','fontsize',14);
    end
    if J > 3
        subplot(2,2,3);Gk = G{J-1};
        Gk =  Gk - spdiags(diag(Gk),0,size(Gk,1),size(Gk,1));
        Gk = tril(Gk);Gk = graph(Gk + Gk');
        plot(Gk,'Layout','force','Iterations',50);axis off;
        title(['Level = ','$J-1$',' = ',num2str(J-1)],'Interpreter','latex','fontsize',14);
        
        subplot(2,2,4);Gk = G{J};
        Gk =  Gk - spdiags(diag(Gk),0,size(Gk,1),size(Gk,1));
        Gk = tril(Gk);Gk = graph(Gk + Gk');
        plot(Gk,'Layout','force','Iterations',50);axis off;
        title(['Level = ','$J$',' = ',num2str(J)],'Interpreter','latex','fontsize',14);
    end
end
end