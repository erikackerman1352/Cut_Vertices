function cutVertices(G)
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

T = dfsSpanningTree(G, 2);

K = [];

if (length(neighbors(T, 1) > 1)
    K(end+1) = [r];
end

for i = 1:numnodes(T)
    non_tree_edge = 
    endpoints = G.Edges.EndNodes(non_tree_edge,:);
    endpoints = findnode(G,{endpoints{1} endpoints{2}});        
    x = min(endpoints); 
    m_w = min(T.Nodes.dfN(x);
    low(T.Nodes.origId(i)) = min(T.Nodes.dfN(i), m_w); 
end
end
