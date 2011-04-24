clear;
%% Build the factors
rows = 100;
cols = 100;
states = 5;
lambdaSmooth = 1.5; % Laplace smoothing parameter
noiseP = 0.3; % proportion of randomly sampled values
[factors, img, noisy_img] = ...
   make_grid_model(rows, cols, states, lambdaSmooth, noiseP);

%% Call the mex interface
compile_all; rehash;
disp('Now you can  run experiment');

%%
%options.alg_type = 'CHROMATIC';
options.alg_type = 'SPLASH';
options.nsamples = 2;
options.nskip = 3;
options.tskip = 5;
options.ncpus = 8;
disp('---------------');
[samples, blfs, nsamples, nchanges] = ...
  sampler(factors, options);

%% Take the last set of beliefs and compute the exected value

pred_img = reshape(cellfun(@(x) (1:length(x)) * x, blfs(:,end)), ...
  rows, cols);
figure(1);subplot(1,3,3); colormap('gray');
imagesc(pred_img); title('Expected Pixel Marginals');