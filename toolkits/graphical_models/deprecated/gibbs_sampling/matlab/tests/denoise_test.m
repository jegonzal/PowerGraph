clear;
%% Build the factors
rows = 100;
cols = 100;
states = 5;
lambdaSmooth = 1.5; % Laplace smoothing parameter
noiseP = 0.3; % proportion of randomly sampled values
[factors, img, noisy_img] = ...
  make_grid_model(rows, cols, states, lambdaSmooth, noiseP);

%%
options.alg_type = 'CHROMATIC';
options.nsamples = 50;
options.nskip = 10;
options.ncpus = 4;
disp('---------------');
[samples, nupdates, nchanges, marginals] = ...
  gibbs_sampler(factors, options);

%% Take the last set of beliefs and compute the exected value
for i = 1:options.nsamples
  %   pred_img = reshape(cellfun(@(x) (1:length(x)) * x, marginals(:,i)), ...
  %     rows, cols);
  pred_img = reshape(arrayfun(@(x) (1:length(x)) * x, samples(:,i)), ...
    rows, cols);
  
  figure(1);subplot(1,3,3); colormap('gray');
  imagesc(pred_img); title(['Sample Image ', num2str(i)]);
end

%% Render the final marginal expectations

pred_img = reshape(cellfun(@(x) (1:length(x)) * x, marginals(:,i)), ...
  rows, cols);
figure(2); colormap('gray');
imagesc(pred_img); title('Expected Pixel Marginals');
