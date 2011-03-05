%% Cut random graph

E = sparse(rand(100, 100) > 0.5);
E = E - diag(diag(E));
E = E + E';
v = ones(100,1);

cut = metiscut(E,v,5);

%% Cut gibbs graph
m = 200;
verts = 1:(m*m);
u = reshape(verts,m,m);
edges = ...
   [reshape(u(2:end, :), m*(m-1),1),...
   reshape(u(1:(end-1), :), m*(m-1), 1);
   reshape(u(:, 2:end), m*(m-1),1),...
   reshape(u(:, 1:(end-1)), m*(m-1), 1)];
E = sparse(edges(:,1), edges(:,2), ones(length(edges),1), m*m, m*m);
E = E' + E;
v = ones(m*m, 1);

%%
cut = metiscut(E,v,50);
imagesc(reshape(cut,m,m));