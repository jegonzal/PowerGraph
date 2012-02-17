%% Make random undirected graph

nvertices = 10000;
approx_edges = 6*nvertices;
density = approx_edges / (nvertices^2) ;

%% create random edges
E = sprand(nvertices, nvertices, density);

spy(E);

%% save random graph
prefix = 'random_graph';
[u,v] = find(E);
tbl = unique([u-1 v-1; v-1 u-1], 'rows');
dlmwrite([prefix, '.tsv'], tbl, 'delimiter', '\t', 'precision', '%8d');