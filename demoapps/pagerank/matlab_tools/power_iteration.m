%clear;
edges = load('~/Documents/graphlabapi/release/demoapps/pagerank/edges.tsv');
nverts = max(edges(:,1)) + 1;
emat = sparse(edges(:,1)+1 , edges(:,2)+1, edges(:,3), nverts, nverts);



%%
alpha = 0.15;
niters = 33;
pr = ones(nverts,1)/nverts;
emat_tr = emat';


for i = 1:niters
    pr_new = (1-alpha) * emat_tr * pr + alpha / nverts; 
    resid = max(abs(pr_new - pr));
    disp(resid);
    pr = pr_new;
end
