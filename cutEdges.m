function [L] = cutEdges(G)

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



L = [];
H = G;
% when the set of frontier edges is not empty
while (numnodes(H) > 1)
    
   index = find(strcmp(H.Nodes.Name, H.Nodes.Name(1)));
   [t, path] = dfsPath(H, index);

   if (degree(H,t) ==1)
       cut_edge = outedges(H, t);
       cut_edge = str2num(cell2mat(H.Edges.Name(cut_edge)));
       L(end+1) = [cut_edge];
       H = rmnode(H, t);                                    
   else
        H = contract(H, t);

   end
   
   
end
L = unique(L);
L = L';
    
end 