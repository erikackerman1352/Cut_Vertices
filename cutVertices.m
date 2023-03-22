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

T = dfsSpanningTree(G, 1);


K = [];

if length(neighbors(T, 1)) >1
    K(end+1) = [1];
end

set1 = [];
set2 = [];
low = [];
for i = 2:numnodes(T)
    s1 = outedges(G, i);
    set1 = cat(2, set1, s1');
    [tf, idx] = ismember(i, T.Nodes.origId);
    s2 = outedges(T, idx);
    s2 = T.Edges.origId(s2);
    set2 = cat(2, set2, s2');
end
    nte = setdiff(set1, set2);

newnte = [];
for j = 1:length(nte)
    endpoints = G.Edges.EndNodes(nte(j),:);
        endpoints = findnode(G,{endpoints{1} endpoints{2}});
        [tf, idx_1] = ismember(endpoints(1), T.Nodes.origId);
        [tf, idx_2] = ismember(endpoints(2), T.Nodes.origId);
        if T.Nodes.dfN(idx_1) == T.Nodes.dfN(idx_2)
            e = nte(j);
            newnte(end+1) = [e];
        end
end


nte = setdiff(nte, newnte);

if (length(neighbors(T,1))==1)
    for i = 2:numnodes(T) 
    
    dfNendp = [];


    for j = 1:length(nte)
        endpoints = G.Edges.EndNodes(nte(j),:);
        endpoints = findnode(G,{endpoints{1} endpoints{2}});
        [tf, idx_1] = ismember(endpoints(1), T.Nodes.origId);
        [tf, idx_2] = ismember(endpoints(2), T.Nodes.origId);
        
        if (T.Nodes.dfN(idx_1) < T.Nodes.dfN(i)) && (T.Nodes.dfN(idx_2) < T.Nodes.dfN(i))
            e = nte(j);
           
        end
         if (T.Nodes.dfN(idx_1) > T.Nodes.dfN(i)) && (T.Nodes.dfN(idx_2) > T.Nodes.dfN(i))
            e = nte(j);
%             nte = setdiff(nte, e);
        end
        if (T.Nodes.dfN(idx_1) == T.Nodes.dfN(i)) || (T.Nodes.dfN(idx_2) == T.Nodes.dfN(i))
            if (T.Nodes.dfN(idx_1) == T.Nodes.dfN(idx_2) )
                dfNendp = [T.Nodes.dfN(i)];
            elseif (T.Nodes.dfN(idx_1) == T.Nodes.dfN(i))
                a = T.Nodes.dfN(idx_1);
                if (T.Nodes.dfN(idx_2) < T.Nodes.dfN(idx_1))
                    dfNendp(end+1) = [T.Nodes.dfN(idx_2)];
                elseif (T.Nodes.dfN(idx_2) > T.Nodes.dfN(idx_1))
                    dfNendp(end+1)  = [T.Nodes.dfN(idx_1)];
                end
            elseif (T.Nodes.dfN(idx_2) == T.Nodes.dfN(i))
                a = T.Nodes.dfN(idx_2);
                if (T.Nodes.dfN(idx_2) < T.Nodes.dfN(idx_1))
                    dfNendp(end+1)  = [T.Nodes.dfN(idx_2)];
                elseif (T.Nodes.dfN(idx_2) > T.Nodes.dfN(idx_1))
                    dfNendp(end+1)  = [T.Nodes.dfN(idx_1)];
                end
            end
            
        end
        if (T.Nodes.dfN(idx_1) > T.Nodes.dfN(i)) && (T.Nodes.dfN(idx_2) < T.Nodes.dfN(i))
            dfNendp(end+ 1) = [T.Nodes.dfN(idx_2)];
        end
        if (T.Nodes.dfN(idx_1) < T.Nodes.dfN(i)) && (T.Nodes.dfN(idx_2) > T.Nodes.dfN(i))
            dfNendp(end+1) = [T.Nodes.dfN(idx_1)];
        end

    end
            low(i) = min(dfNendp);     

    end
    for i = 2:numnodes(T)
    if length(neighbors(T, i)) ==1
        if length(outedges(G, T.Nodes.origId(i))) <2
            low(i) = T.Nodes.dfN(i);
        end
    end
    end
    
end

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
% if length(neighbors(T, 1)) == 1
% for i = 2:(numnodes(T)-1)
%     if (low(i+1) >= T.Nodes.dfN(i))
%         K(end+1) = [T.Nodes.origId(i)];
%     end
% end
% K = unique(K);
% K = K';
% end
if isempty(K)
    K = [];
end
end
