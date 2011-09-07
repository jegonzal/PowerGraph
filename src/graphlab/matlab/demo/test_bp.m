%%
% This is an image denoising example with BP
% We construct a noisy image of a rainbow and connects the
% pixels in a grid pattern to form a pairwise Markov Random Field.

%% SET THIS!!
% this should be the location of graphlab.a etc
BINARY_DIRECTORY = [getenv('HOME') '/graphlabapi/graphlabapi/release/src/graphlab'];
%%
% Image dimensions and arity
arity = 5;
imgdim = 100;
%% generate images
% This generates a clean rainbow image
cleanimg = zeros(imgdim,imgdim);
maxd = [imgdim,imgdim]/2;
maxd = sqrt(sum(maxd.*maxd));
for i = 1:(imgdim/2)
    for j = 1:imgdim
        d = [i,j] - ([imgdim,imgdim]/2);
        d = sqrt(sum(d.*d));
        cleanimg(i,j) = round(d / maxd * arity);
    end
end
cleanimg(cleanimg >= arity) = arity - 1;

% make noisy image by adding random gaussian noise
noisyimg = cleanimg + randn(imgdim);
%% display images
figure;
image(cleanimg / arity * 100);
colormap('gray');
figure;
image(noisyimg / arity * 100);
colormap('gray');
%% generate vertex data
vdata = {};
for i = 1:imgdim
    for j = 1:imgdim
        vtx.belief = ones(1,arity); % Reset the vertex belief
        % set the unary potential to be a gaussian pdf around the noisy image
        vtx.unary = pdf('norm',0:(arity-1), noisyimg(i,j), 1); 
        % normalize
        vtx.unary = vtx.unary ./ sum(vtx.unary);
        % Vertices are numbered from 1. Assign the vertex data.
        vdata{(i-1)*100+j} = vtx;
    end
end

%% generate laplace edgedata

edge.binary = zeros(arity);
for i = 1:arity
    for j = 1:arity
        edge.binary(i,j) = exp(-3 * abs(i-j));
    end
end
edge.msg = ones(1,5);
% We only need one entry in the edge data since all edges have identical values
edata = {edge};
%% generate adjacency matrix
% connect up data in a grid
adj =sparse(length(vdata),length(vdata));
for xi = 1:imgdim
    for xj = 1:imgdim
        if (xi - 1 >= 1) 
            adj((xi-1)*100+xj, (xi-2)*100+xj) = 1;
        end
        if (xi + 1 <= imgdim) 
            adj((xi-1)*100+xj, (xi)*100+xj) = 1;
        end
        if (xj - 1 >= 1) 
            adj((xi-1)*100+xj, (xi-1)*100+xj-1) = 1;
        end
        if (xj + 1 <= imgdim) 
            adj((xi-1)*100+xj, (xi-1)*100+xj+1) = 1;
        end
    end
end
%% compile update function
compile_update_function({'bp_update'}, vdata{1},edata{1}, BINARY_DIRECTORY, 'bp', 'bp', 3);
%% set options
options.initial_schedule(1).update_function = 'bp_update';
options.initial_schedule(1).vertices=uint32(1:(imgdim * imgdim));
options.initial_schedule(1).priorities=ones(size(options.initial_schedule(1).vertices));
options.scheduler = 'splash(50)';
options.ncpus = 4;
options.scope = 'edge';
%% cd in the BP directory and call
cd bp
%%
[v2,adj2,e2] = bp(vdata,adj,edata, options);
%%
cd ..
%% display output

outputimg = zeros(imgdim);
for i = 1:imgdim
    for j = 1:imgdim
        [c,idx] = max(v2{(i-1)*100+j}.belief);
        outputimg(i,j) = idx;
    end
end
figure;
imagesc(outputimg);
colormap('gray');
