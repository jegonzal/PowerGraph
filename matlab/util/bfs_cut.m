function asgs = bfs_cut(E, block_count) 

% Make E symmetric
E = E' + E;

unassigned = 1:length(E);
assigned = [];

asgs = -1*ones(length(unassigned),1);

maxblocksize = (length(unassigned) / block_count) + 1;

progress = 0;
owner = 1;


while(not(isempty(unassigned))) 
   asg_count = 0;
   visited = [];
   queue = []; 
   while(asg_count < maxblocksize && not(isempty(unassigned))) 
      if(isempty(queue))
         queue = unassigned(1);
         visited = queue;
      end
      
      vertex = queue(1); queue(1) = [];
      asgs(vertex) = owner;
      assigned = [assigned, vertex];
      unassigned = setdiff(unassigned, vertex);
      asg_count = asg_count + 1;
      
      current_progress = round(100*(length(asgs) - length(unassigned)) / length(asgs));
      if(current_progress ~= progress)
         progress = current_progress;
         disp(['Progress: ', num2str(progress)]);
      end
      
      
      % add all neighbors
      neighbors = find(E(vertex,:));
      valid_neighbors = setdiff(setdiff(neighbors, assigned), visited);
      queue = [queue, valid_neighbors];
      visited = union(queue, valid_neighbors);
   end
   
   
   
   disp(['Owner : ', num2str(owner)]);
   owner = owner + 1;
end




end







