%clear;
edges = load('~/Documents/graphlabapi/release/demoapps/pagerank/edges.tsv');
%%
nverts = max(max(edges(:,1:2))) + 1;
emat = sparse(edges(:,1)+1 , edges(:,2)+1, 1, nverts, nverts);

%% NOrmalize
Z  = full(sum(emat,2));
emat = sparse(1:nverts, 1:nverts, 1./Z, nverts,nverts) * emat;


%%
alpha = 0.15;
pr = ones(nverts,1)/nverts;
emat_tr = emat';

%%
%pr_new = (1-alpha) * emat_tr * pr + alpha / nverts;
pr_new = (1-alpha) * emat * pr + alpha / nverts;
sum(pr_new)
resid = max(abs(pr_new - pr));
disp(resid);
pr = pr_new;

%%
niters = 33;
for i = 1:niters
    pr_new = (1-alpha) * emat_tr * pr + alpha / nverts; 
    resid = max(abs(pr_new - pr));
    disp(resid);
    pr = pr_new;
end
