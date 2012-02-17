nvertices = 10000;
approx_edges = 3*nvertices;
density = approx_edges / (nvertices^2) ;

%% create random edges
E = sprand(nvertices, nvertices, density);

spy(E);

%% save random graph
prefix = 'random_graph';
[u,v,w] = find(E);
tbl = [u-1,v-1];
dlmwrite([prefix, '.tsv'], tbl, 'delimiter', '\t', 'precision', '%8d');