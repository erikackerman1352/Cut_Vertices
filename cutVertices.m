function [K] = cutVertices(G)
% check if vertices have names
if (~sum(ismember(G.Nodes.Properties.VariableNames,'Name')))
    % if not, give names using its indices
    Vnames = int2str(1:numnodes(G));
    G.Nodes.Name = split(Vnames);
end

% check if edges have names
if (~sum(ismember(G.Edges.Properties.VariableNames,'Name')))
    % if not, give names using its indices
    Enames = int2str(1:numedges(G));
    G.Edges.Name = split(Enames);
end

% generate dfs spanning tree from the 1st vertex
T = dfsSpanningTree(G, 1);

% set an empty set K
K = [];

% if the children of T is more than 1, then add vertex 1 to K
if length(neighbors(T, 1)) >1
    K(end+1) = [1];
end

%the next 2 for loops will find the edges of in T (without the non root edges from G)
set1 = [];
set2 = [];
low = [];

% loops through all nodes of T
for i = 2:numnodes(T)
    
    % find the outedges at each nodes, but in G
    s1 = outedges(G, i);
    set1 = cat(2, set1, s1');
    
    % find the index of each nodes in T 
    [tf, idx] = ismember(i, T.Nodes.origId);
    
    % find the outedges at each nodes, but in T
    s2 = outedges(T, idx);
    s2 = T.Edges.origId(s2);
    set2 = cat(2, set2, s2');
end

   % the nte set is the difference between edges in T and G, which are none root edges
    nte = setdiff(set1, set2);

newnte = [];

% loop through all edges in nte
for j = 1:length(nte)
        
    %find the endpoints of those edges
        endpoints = G.Edges.EndNodes(nte(j),:);
        endpoints = findnode(G,{endpoints{1} endpoints{2}});
        
        %find the indexes of those points 
        [tf, idx_1] = ismember(endpoints(1), T.Nodes.origId);
        [tf, idx_2] = ismember(endpoints(2), T.Nodes.origId);
        
        %if the dfN at both endpoints are the same, this is a self loop edge.
        if T.Nodes.dfN(idx_1) == T.Nodes.dfN(idx_2)
            e = nte(j);
            newnte(end+1) = [e];
        end
end

% eliminate those self loops from nte
nte = setdiff(nte, newnte);


% For this problem, I seperate the process in 2 cases,
% 1 case with node 1 has 1 children, and 1 case with node 1 has more than 1 children

% 1st case, T has 1 children
if (length(neighbors(T,1))==1)
    
    %loop throgh all nodes of T
    for i = 2:numnodes(T) 
        dfNendp = [];

    % for each nodes in T, analyze all non-tree edges in nte.  
    % the idea is analyze the dfN of the endpoint.
    
    % the process is seperated in many cases. 
    
    % case 1: if both of the dfN of the endpoints of such edge is SMALLER than the dfN of the node 
    % we are analyzing, say node i, this means the edge is connecting the nodes ABOVE the node i, 
    % which is has no affect on determining the low number => do nothing
    
    % case 2: if both of the dfN of the endpoints of such edge is LARGER than the dfN of the node 
    % we are analyzing, say node i, this means the edge is connecting the nodes BELOW the node i, 
    % which is has no affect on determining the low number => do nothing
    
    % case 3: if either dfN of one of the endpoint equal the dfN of the node i, this will effect the low number of node i
    % This effect the low number
    
        % seperate into smaller senarios:
        % senario 1: if one of two endpoints equal the dfN of the node i, then the nontree root is starting from node i
        % if the node of the other endpoint is GREATER than dfN of node i, then that node is BELOW node i in tree T, I 
        % record the dfN of node i
        
        % senario 2: if one of two endpoints equal the dfN of the node i, then the nontree root is starting from node i
        % if the node of the other endpoint is SMALLER than dfN of node i, then that node is ABOVE node i in tree T, I 
        % record the dfN of the other node (not i).
        
   % case 4: if dfN one of the endpoint is greater and the other is smaller than dfN of node i, then we choose dfN of 
   % the smaller one to be saved.
   
    % After looping through all the nodes, I save them into dfNendp, and pick the minimum one. It is the low number of that node i
        
    for j = 1:length(nte)
        
        %find endpoints of those edges
        endpoints = G.Edges.EndNodes(nte(j),:);
        endpoints = findnode(G,{endpoints{1} endpoints{2}});
        [tf, idx_1] = ismember(endpoints(1), T.Nodes.origId);
        [tf, idx_2] = ismember(endpoints(2), T.Nodes.origId);
        
        %case 1
        if (T.Nodes.dfN(idx_1) < T.Nodes.dfN(i)) && (T.Nodes.dfN(idx_2) < T.Nodes.dfN(i))
             
        end
        
        %case2
         if (T.Nodes.dfN(idx_1) > T.Nodes.dfN(i)) && (T.Nodes.dfN(idx_2) > T.Nodes.dfN(i))
           
         end
         
         % case 3
        if (T.Nodes.dfN(idx_1) == T.Nodes.dfN(i)) || (T.Nodes.dfN(idx_2) == T.Nodes.dfN(i))
            if (T.Nodes.dfN(idx_1) == T.Nodes.dfN(idx_2) )
                dfNendp = [T.Nodes.dfN(i)];
                
                %senario 1
            elseif (T.Nodes.dfN(idx_1) == T.Nodes.dfN(i))
                a = T.Nodes.dfN(idx_1);
                if (T.Nodes.dfN(idx_2) < T.Nodes.dfN(idx_1))
                    dfNendp(end+1) = [T.Nodes.dfN(idx_2)];
                elseif (T.Nodes.dfN(idx_2) > T.Nodes.dfN(idx_1))
                    dfNendp(end+1)  = [T.Nodes.dfN(idx_1)];
                end
                
                %senario 2
            elseif (T.Nodes.dfN(idx_2) == T.Nodes.dfN(i))
                a = T.Nodes.dfN(idx_2);
                if (T.Nodes.dfN(idx_2) < T.Nodes.dfN(idx_1))
                    dfNendp(end+1)  = [T.Nodes.dfN(idx_2)];
                elseif (T.Nodes.dfN(idx_2) > T.Nodes.dfN(idx_1))
                    dfNendp(end+1)  = [T.Nodes.dfN(idx_1)];
                end
            end
        end
        
        %case 4
        if (T.Nodes.dfN(idx_1) > T.Nodes.dfN(i)) && (T.Nodes.dfN(idx_2) < T.Nodes.dfN(i))
            dfNendp(end+ 1) = [T.Nodes.dfN(idx_2)];
        end
        if (T.Nodes.dfN(idx_1) < T.Nodes.dfN(i)) && (T.Nodes.dfN(idx_2) > T.Nodes.dfN(i))
            dfNendp(end+1) = [T.Nodes.dfN(idx_1)];
        end

    end
            low(i) = min(dfNendp);     

    end
    
    %case 5: for the edge that goes to the nodes that are also leaves
    for i = 2:numnodes(T)
    if length(neighbors(T, i)) ==1
        if length(outedges(G, T.Nodes.origId(i))) <2
            low(i) = T.Nodes.dfN(i);
        end
    end
    end
    
end


% case if the node 1 has more than 1 childs, which leads to multiple branched
% I have to work on each each branch the same way as above. However,
% the main difference is isolating the non-tree edges for that belongs on each branch.

if (length(neighbors(T, 1))>1)
edge = [];
ns = neighbors(T,1);
ns_ns1 = neighbors(G, T.Nodes.origId(ns(1)));
newns = [];
for i = 1:length(ns_ns1)
    [tf, idx_1] = ismember(ns_ns1(i), T.Nodes.origId);
    if (T.Nodes.origId(ns_ns1(i) ~=1))
        ns_ns_2 = neighbors(T, idx_1);
        newns = cat(2, newns, ns_ns_2');
        
    end
end
for k = 1:length(newns)
    n = neighbors(T, newns(k));
    newns = cat(2, newns, n');
    newns = unique(newns);
    leftnode = newns(newns~=1);
end

for i = 1:length(leftnode)
    e = outedges(T, leftnode(i));
    edge = cat(2, edge, e');
    edge = unique(edge);
end
edge2 = [];
for i = 1:length(leftnode)
    e = outedges(G, leftnode(i));
    edge2 = cat(2, edge2, e');
    edge2 = unique(edge2);
end
    edge1 = T.Edges.origId(edge)';
    newnte = setdiff(edge2, edge1);
    T_leftnode = [];
    for i = 1:length(leftnode)
        [tf, id] = ismember(leftnode(i), T.Nodes.origId);
        T_leftnode(end+1) = [id];
    end
for i = 1:length(T_leftnode)
    
    dfNendp = [];


    for j = 1:length(newnte)
        endpoints = G.Edges.EndNodes(newnte(j),:);
        endpoints = findnode(G,{endpoints{1} endpoints{2}});
        [tf, idx_1] = ismember(endpoints(1), T.Nodes.origId);
        [tf, idx_2] = ismember(endpoints(2), T.Nodes.origId);
        
        if (T.Nodes.dfN(idx_1) < T.Nodes.dfN(T_leftnode(i))) && (T.Nodes.dfN(idx_2) < T.Nodes.dfN(T_leftnode(i)))
            e = newnte(j);
           
        end
         if (T.Nodes.dfN(idx_1) > T.Nodes.dfN(T_leftnode(i))) && (T.Nodes.dfN(idx_2) > T.Nodes.dfN(T_leftnode(i)))
            e = newnte(j);
%             nte = setdiff(nte, e);
        end
        if (T.Nodes.dfN(idx_1) == T.Nodes.dfN(T_leftnode(i))) || (T.Nodes.dfN(idx_2) == T.Nodes.dfN(T_leftnode(i)))
            if (T.Nodes.dfN(idx_1) == T.Nodes.dfN(idx_2) )
                dfNendp = [T.Nodes.dfN(T_leftnode(i))];
            elseif (T.Nodes.dfN(idx_1) == T.Nodes.dfN(T_leftnode(i)))
                a = T.Nodes.dfN(idx_1);
                if (T.Nodes.dfN(idx_2) < T.Nodes.dfN(idx_1))
                    dfNendp(end+1) = [T.Nodes.dfN(idx_2)];
                elseif (T.Nodes.dfN(idx_2) > T.Nodes.dfN(idx_1))
                    dfNendp(end+1)  = [T.Nodes.dfN(idx_1)];
                end
            elseif (T.Nodes.dfN(idx_2) == T.Nodes.dfN(T_leftnode(i)))
                a = T.Nodes.dfN(idx_2);
                if (T.Nodes.dfN(idx_2) < T.Nodes.dfN(idx_1))
                    dfNendp(end+1)  = [T.Nodes.dfN(idx_2)];
                elseif (T.Nodes.dfN(idx_2) > T.Nodes.dfN(idx_1))
                    dfNendp(end+1)  = [T.Nodes.dfN(idx_1)];
                end
            end
            
        end
        if (T.Nodes.dfN(idx_1) > T.Nodes.dfN(T_leftnode(i))) && (T.Nodes.dfN(idx_2) < T.Nodes.dfN(T_leftnode(i)))
            dfNendp(end+ 1) = [T.Nodes.dfN(idx_2)];
        end
        if (T.Nodes.dfN(idx_1) < T.Nodes.dfN(T_leftnode(i))) && (T.Nodes.dfN(idx_2) > T.Nodes.dfN(T_leftnode(i)))
            dfNendp(end+1) = [T.Nodes.dfN(idx_1)];
        end

    end
    if isempty(dfNendp)
        low(T_leftnode(i)) = 0;
    else
        low(T_leftnode(i)) = min(dfNendp);

    end
    for i = 2:numnodes(T)
    if length(neighbors(T, i)) ==1
        if length(outedges(G, T.Nodes.origId(i))) <2
            low(i) = T.Nodes.dfN(i);
        end
    end
end
end
    
edge = [];
ns = neighbors(T,1);
ns_ns1 = neighbors(G, T.Nodes.origId(ns(2)));
newns = [];
for i = 1:length(ns_ns1)
    [tf, idx_1] = ismember(ns_ns1(i), T.Nodes.origId);
    if (T.Nodes.origId(ns_ns1(i) ~=1))
        ns_ns_2 = neighbors(T, idx_1);
        newns = cat(2, newns, ns_ns_2');
        
    end
end
for k = 1:length(newns)
    n = neighbors(T, newns(k));
    newns = cat(2, newns, n');
    newns = unique(newns);
    rightnode = newns(newns~=1);
end

for i = 1:length(rightnode)
    e = outedges(T, rightnode(i));
    edge = cat(2, edge, e');
    edge = unique(edge);
end
edge2 = [];
for i = 1:length(rightnode)
    e = outedges(G, rightnode(i));
    edge2 = cat(2, edge2, e');
    edge2 = unique(edge2);
end
    edge1 = T.Edges.origId(edge)';
    newnte = setdiff(edge2, edge1);
    T_rightnode = [];
    for i = 1:length(rightnode)
        [tf, id] = ismember(rightnode(i), T.Nodes.origId);
        T_rightnode(end+1) = [id];
    end
    for i = 1:length(T_rightnode)
    
    dfNendp = [];


    for j = 1:length(newnte)
        endpoints = G.Edges.EndNodes(newnte(j),:);
        endpoints = findnode(G,{endpoints{1} endpoints{2}});
        [tf, idx_1] = ismember(endpoints(1), T.Nodes.origId);
        [tf, idx_2] = ismember(endpoints(2), T.Nodes.origId);
        
        if (T.Nodes.dfN(idx_1) < T.Nodes.dfN(T_rightnode(i))) && (T.Nodes.dfN(idx_2) < T.Nodes.dfN(T_rightnode(i)))
            e = newnte(j);
           
        end
         if (T.Nodes.dfN(idx_1) > T.Nodes.dfN(T_rightnode(i))) && (T.Nodes.dfN(idx_2) > T.Nodes.dfN(T_rightnode(i)))
            e = newnte(j);
%             nte = setdiff(nte, e);
        end
        if (T.Nodes.dfN(idx_1) == T.Nodes.dfN(T_rightnode(i))) || (T.Nodes.dfN(idx_2) == T.Nodes.dfN(T_rightnode(i)))
            if (T.Nodes.dfN(idx_1) == T.Nodes.dfN(idx_2) )
                dfNendp = [T.Nodes.dfN(T_rightnode(i))];
            elseif (T.Nodes.dfN(idx_1) == T.Nodes.dfN(T_rightnode(i)))
                a = T.Nodes.dfN(idx_1);
                if (T.Nodes.dfN(idx_2) < T.Nodes.dfN(idx_1))
                    dfNendp(end+1) = [T.Nodes.dfN(idx_2)];
                elseif (T.Nodes.dfN(idx_2) > T.Nodes.dfN(idx_1))
                    dfNendp(end+1)  = [T.Nodes.dfN(idx_1)];
                end
            elseif (T.Nodes.dfN(idx_2) == T.Nodes.dfN(T_rightnode(i)))
                a = T.Nodes.dfN(idx_2);
                if (T.Nodes.dfN(idx_2) < T.Nodes.dfN(idx_1))
                    dfNendp(end+1)  = [T.Nodes.dfN(idx_2)];
                elseif (T.Nodes.dfN(idx_2) > T.Nodes.dfN(idx_1))
                    dfNendp(end+1)  = [T.Nodes.dfN(idx_1)];
                end
            end
            
        end
        if (T.Nodes.dfN(idx_1) > T.Nodes.dfN(T_rightnode(i))) && (T.Nodes.dfN(idx_2) < T.Nodes.dfN(T_rightnode(i)))
            dfNendp(end+ 1) = [T.Nodes.dfN(idx_2)];
        end
        if (T.Nodes.dfN(idx_1) < T.Nodes.dfN(T_rightnode(i))) && (T.Nodes.dfN(idx_2) > T.Nodes.dfN(T_rightnode(i)))
            dfNendp(end+1) = [T.Nodes.dfN(idx_1)];
        end

    end
            if isempty(dfNendp)
        low(T_rightnode(i)) = 0;
    else
        low(T_rightnode(i)) = min(dfNendp);

    end
    end
end

if length(neighbors(T, 1)) >1 || length(neighbors(T,1)) ==1
for i = 2:(numnodes(T)-1)
    e = outedges(T, i);
    for j = 1:length(e)
        idx = T.Edges.EndNodes(e(j),:);
        idx = findnode(G,{idx{1} idx{2}});
        [tf, endpoints1] = ismember(idx(1), T.Nodes.origId);
        [tf, endpoints2] = ismember(idx(2), T.Nodes.origId);
        if T.Nodes.dfN(endpoints1) == T.Nodes.dfN(i)
            if T.Nodes.dfN(endpoints2) < T.Nodes.dfN(i)
            elseif T.Nodes.dfN(endpoints2) > T.Nodes.dfN(i)
                if(low(endpoints2) >= T.Nodes.dfN(i))
                    K(end+1) = [T.Nodes.origId(i)];
                end
            end
        elseif T.Nodes.dfN(endpoints2) == T.Nodes.dfN(i)
            if T.Nodes.dfN(endpoints1) < T.Nodes.dfN(i)
            elseif T.Nodes.dfN(endpoints1) > T.Nodes.dfN(i)
                if(low(endpoints1) >= T.Nodes.dfN(i))
                    K(end+1) = [T.Nodes.origId(i)];
                end
            end
        elseif T.Nodes.dfN(endpoints1) < T.Nodes.dfN(i) || T.Nodes.dfN(endpoints2) < T.Nodes.dfN(i)
        elseif T.Nodes.dfN(endpoints1) > T.Nodes.dfN(i) || T.Nodes.dfN(endpoints2) > T.Nodes.dfN(i)
            if T.Nodes.dfN(endpoints2) > T.Nodes.dfN(i)
                if(low(endpoints2) >= T.Nodes.dfN(i))
                    K(end+1) = [T.Nodes.origId(i)];
                end
            elseif T.Nodes.dfN(endpoints1) > T.Nodes.dfN(i)
                if(low(endpoints1) >= T.Nodes.dfN(i))
                    K(end+1) = [T.Nodes.origId(i)];
                end
            end
        end
    end
end

K = unique(K);
K = K';
        
end

if isempty(K)
    K = [];
end
end % end function cutVertices