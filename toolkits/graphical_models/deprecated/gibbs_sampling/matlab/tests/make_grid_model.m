%% This code generates a grid model
function [factors, img, noisy_img] = make_grid_model(rows, cols, states, ...
   lambdaSmooth, noiseP)

% Create a virtual image
[u,v] = meshgrid(linspace(0,1,rows), linspace(0,1,cols));
img = (1 + cos(1./sqrt((u-.5).^2 + (v-.5).^2)) )/2 + u.^2;
img = (img - min(img(:)))/(max(img(:)) - min(img(:)));
img = (states - 1) * img;
img = round(img) + 1;
figure(1); clf();subplot(1,3,1); colormap('gray');
imagesc(img);
title('Original Image');

% add noise
mask = rand(rows,cols) < noiseP;
noise = ceil(states*rand(rows,cols));
noisy_img = mask .* noise + ~mask .* img;
figure(1); subplot(1,3,2); colormap('gray');
imagesc(noisy_img);
title('Noisy Image');
% Build the edge factor table (in log form)
[u,v] = meshgrid(1:states, 1:states);
edgetbl = -lambdaSmooth * abs(u - v);


% Build the node factor tableS based on the noise model
nodetbls = zeros(rows*cols, states) + ...
   noiseP/(states - 1);
ind = sub2ind([rows * cols, states], (1:(rows*cols))', noisy_img(:));
nodetbls(ind) = 1-noiseP;
nodetbls = log(nodetbls);


% Get all the edges and variables
vars = 1:(rows*cols);
gridvars = reshape(vars, rows, cols);
edges = ...
   [reshape(gridvars(1:(end-1),:), (rows-1) * cols,1), ...
   reshape(gridvars(2:end,:), (rows-1) * cols,1); ...
   reshape(gridvars(:,1:(end-1)), rows * (cols-1), 1), ...
   reshape(gridvars(:,2:end), rows * (cols-1),1)];


% construct the actual factors
factors = cell(length(vars) + length(edges), 1);
index = 1;
for i = 1:length(vars)
   factors{index} = table_factor(vars(i), nodetbls(vars(i),:));   
   index = index + 1;
end
index = length(vars)+1;
for i = 1:length(edges)
   factors{index} = table_factor(sort(edges(i,:)), edgetbl);   
   index = index + 1;
end



end