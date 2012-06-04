clear;
RandStream.setDefaultStream(RandStream.create('mt19937ar', 'seed', 5849))

smoothing = 0.000001;

var_sizes = [3,5,3,2,4];
vars = [1,2,3];
factors{1} = table_factor(vars, log(rand(var_sizes(vars)) + smoothing));

vars = [2,3,4];
factors{2} = table_factor(vars, log(rand(var_sizes(vars)) + smoothing));

vars = [1,3,5];
factors{3} = table_factor(vars, log(rand(var_sizes(vars)) + smoothing));

vars = [4,5];
factors{4} = table_factor(vars, log(rand(var_sizes(vars)) + smoothing));


maxasg = prod(var_sizes);

joint = zeros(var_sizes);

%% compute joint
for i = 1:maxasg
  asg = ind2asg(var_sizes, i);
  for j = 1:length(factors)
    subi = asg2ind(var_sizes(factors{j}.vars), asg(factors{j}.vars));
    joint(i) = joint(i)  + factors{j}.logP(subi);
  end
end

P = exp(joint) / sum(exp(joint(:)));


%% run the Chromatic Sampler
options.alg_type = 'CHROMATIC';
options.nsamples = 1000;
options.nskip = 100;
options.ncpus = 1;

samples = gibbs_sampler(factors, options);

P_est = zeros(var_sizes);
for i = 1:options.nsamples
  ind = asg2ind(var_sizes, samples(:,i)');
  P_est(ind) = P_est(ind) + 1;
end
P_est = P_est ./ sum(P_est(:));


error = abs(P_est - P);
disp(['Chromatic error: ', num2str(max(error(:)))]);

%% run the sampler
options.alg_type = 'SPLASH';
options.nsamples = 1000;
options.nskip = 100;
options.ncpus = 1;
options.treewidth=5;
samples = gibbs_sampler(factors, options);


P_est = zeros(var_sizes);
for i = 1:options.nsamples
  ind = asg2ind(var_sizes, samples(:,i)');
  P_est(ind) = P_est(ind) + 1;
end
P_est = P_est ./ sum(P_est(:));


error = abs(P_est - P);
disp(['Splash error: ', num2str(max(error(:)))]);