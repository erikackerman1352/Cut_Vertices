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
for i = 2:numnodes(T) 
    [tf, idx] = ismember(i, T.Nodes.origId);
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
for i = 2:(numnodes(T)-1)
    if (low(i+1) >= T.Nodes.dfN(i))
        K(end+1) = [T.Nodes.origId(i)];
    end
end

   
   
   
end
