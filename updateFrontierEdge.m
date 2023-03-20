function S_new = updateFrontierEdge(G,S)

%First, sort S from small to big
sort_S = sort(S);

S_new = [];
%Start a loop across the length of Nodes of G.
  for i = 1:length(G.Nodes.Name)
 
      % Look for the node with the highest dfN, which is also the newest
      % discovered node
      if (G.Nodes.dfN(i) == max(G.Nodes.dfN))
      
        % Set S_new as the frontier edges from that node. However, this
        % S includes both discovred and undiscovered edges
        [S_new,nV] = outedges(G,G.Nodes.Name(i));
%         S_new(end+1) = cat(2, S_new, newS');
        

        end
  end

%Now filter out the discovered edges
%set an empty discovered S for later
discoveredS = [];

%start a loop across the length of the S_new found above
for k = 1:length(S_new)
        
        %locate the endpoints of each edge at each k-th
        endpoints = G.Edges.EndNodes(S_new(k),:);
        endpoints = findnode(G,{endpoints{1} endpoints{2}});
        
        %if the dfN of both of the endpoints of that k-th edge is greater
        %or equal to 0, then that edges is discovred.
        if (G.Nodes.dfN(endpoints(1))>=0) && (G.Nodes.dfN(endpoints(2))>=0)
            
            %call the k-th discovered edge is newS
            newS = (S_new(k));
            
            %append all the discovered edges to discoveredS
            discoveredS(end+1) = [newS];
    
        end
        

end

%use setdiff to extract the undiscovered from S_new vs the discovredS.
S_new = setdiff(S_new, discoveredS);

%From this line down, it is used to check for any nodes that have not been discovered
%set an empty end points, endp, for later use
endp = [];

% if S is empty, but in order to be sure all nodes are discovered
if (isempty(S_new))
    
    %run a loop across all nodes again
    for j = 1:length(G.Nodes.Name)
        
        %if there exist any nodes that has dfN is Inf
        if (isinf(G.Nodes.dfN(j)))
            
            % extract all edges from that nodes
            S_new = outedges(G, j);
            
            %looping through all the edges in that S_new
            for k = 1:length(S_new)
                
                %extract the endpoints of those edges
                endpoints = G.Edges.EndNodes(S_new(k),:);
                endpoints = findnode(G,{endpoints{1} endpoints{2}});
                
                %append those nodes into 1 list, and only keep the unique elements of them
                endp = cat(2, endp, endpoints');
                endp = unique(endp);
            end
            
           
            
        
    
    e_all = [];
    for i = 1:length(endp)
        ns = neighbors(G, endp(i));
        for k = 1:length(ns)
            e = findedge(G, endp(i), ns(k));
            if length(e)>=1;
                for j  =  1:length(e)
                    e_all(end+1) = [e(j)];
                end
            else
            e_all(end+1) = [e];
            end
            
        end
        e_all = unique(e_all);
        e_all = e_all(e_all~=0);
    end
    discoveredS = [];

%start a loop across the length of the S_new found above
for k = 1:length(e_all)
        
        %locate the endpoints of each edge at each k-th
        endpoints = G.Edges.EndNodes(e_all(k),:);
        endpoints = findnode(G,{endpoints{1} endpoints{2}});
        
        %if the dfN of both of the endpoints of that k-th edge is greater
        %or equal to 0, then that edges is discovred.
        if (G.Nodes.dfN(endpoints(1))>=0) && (G.Nodes.dfN(endpoints(2))>=0)
            
            %call the k-th discovered edge is newS
            newS = (e_all(k));
            
            %append all the discovered edges to discoveredS
            discoveredS(end+1) = [newS];
    
        end
        

end

%use setdiff to extract the undiscovered from S_new vs the discovredS.
S_new = setdiff(e_all, discoveredS);
undiscoveredS = [];
for k = 1:length(S_new)
     %locate the endpoints of each edge at each k-th
        endpoints = G.Edges.EndNodes(S_new(k),:);
        endpoints = findnode(G,{endpoints{1} endpoints{2}});
        if (G.Nodes.dfN(endpoints(1))==-Inf) && (G.Nodes.dfN(endpoints(2))==-Inf)
            
            %call the k-th discovered edge is newS
            newS = (S_new(k));
            
            %append all the discovered edges to discoveredS
            undiscoveredS(end+1) = [newS];
    
        end
end


        
S_new = setdiff(S_new, undiscoveredS);
endp = [];
for k = 1:length(S_new)
                
                %extract the endpoints of those edges
                endpoints = G.Edges.EndNodes(S_new(k),:);
                endpoints = findnode(G,{endpoints{1} endpoints{2}});
                
                %append those nodes into 1 list, and only keep the unique elements of them
                endp = cat(2, endp, endpoints');
                endp = unique(endp);
            end
            
            % Find the nodes with the max dfN among those nodes
            [pre_dfN, pos] = max(G.Nodes.dfN(endp));
            
            %Find edge from the node above and the node with the Inf dfN. That
            newS = outedges(G, endp(pos)); 
            S_new = setdiff(newS, discoveredS);

        end
    end
else
        S_new = S_new;

end

end % end function updateFrontierEdge