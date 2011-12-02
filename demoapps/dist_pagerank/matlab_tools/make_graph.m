dir = '/mnt/bigbrofs/usr5/graphlab/testdata/pagerank';
file = 'web-Google.txt';
prefix = 'google';
data = importdata([dir, '/', file]);


sizeof.int = 'uint32'; 
sizeof.real = 'double';

%% Process the data
edges = data.data;
vertices = unique(data.data(:));

%% Rescale
edges = edges + 1;
vertices = vertices + 1;
n = max(vertices);
m = length(edges);

%%  Duplicate edge addressing
[edges, junk, revmap] = unique(edges, 'rows');
eweights = histc( revmap, 1:length(edges) );

%% Make sparse matrix
E = sparse(edges(:,1), edges(:,2), eweights, n, n);
% Add self links
E = E + sparse(1:n, 1:n, ones(n,1), n, n);
% E_ij i -> j count

%% Normalize rows
% E = sparse(1:n , 1:n, 1./sum(E,2), n, n) * E;


%% Save the graph
[u,v,w] = find(E);
tbl = [u-1,v-1,w];
dlmwrite([prefix, '.tsv'], tbl, 'delimiter', '\t', 'precision', '%7d');
