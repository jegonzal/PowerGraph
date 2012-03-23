function evaluate_powerlaw(filename)
%%
data = importdata(filename) + 1;
nverts = max(data(:));
adj = sparse(data(:,1), data(:,2), ones(1,size(data,1)), nverts, nverts);
adj = adj | adj';
figure(1); clf(); spy(adj);

degree = full(sum(adj));
counts = histc(degree, 1:max(degree));
figure(2); clf(); loglog(counts, '-x');
xlabel('log(degree)');
ylabel('log(count)');

% ordered_degree = sort(degree, 'descend');
% figure(2); clf(); loglog(1:nverts, ordered_degree, '-x');


end


