nvertices = 10;
approx_edges = nvertices;
density = approx_edges / (nvertices^2) ;

%% create random edges
E = sprand(nvertices, nvertices, density);
% Add self links
E = E + sparse(1:nvertices, 1:nvertices, ones(nvertices,1), ...
    nvertices, nvertices);

%% Make it symmetric (not strictly necessary)
E = E + E';

%% Normalize rows
E = sparse(1:nvertices , 1:nvertices, 1./sum(E,2), ...
    nvertices, nvertices) * E;


%% save random graph
prefix = 'random_graph';
[u,v,w] = find(E);
tbl = [u-1,v-1,w];
dlmwrite([prefix, '.tsv'], tbl, 'delimiter', '\t', 'precision', '%8d');



