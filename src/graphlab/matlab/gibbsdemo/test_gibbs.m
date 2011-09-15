setup_stuff
%% display images
figure;
image(cleanimg / arity * 100);
colormap('gray');
figure;
image(noisyimg / arity * 100);
colormap('gray');


%% compile update function
compile_update_function({'gibbs_update'}, vdata{1},edata{1}, BINARY_DIRECTORY, 'gibbs', 'gibbs', 3);
%% set options
options.initial_schedule(1).update_function = 'gibbs_update';
options.initial_schedule(1).vertices=uint32(1:(imgdim * imgdim));
options.initial_schedule(1).priorities=ones(size(options.initial_schedule(1).vertices));
options.scheduler = 'chromatic(max_iterations=100)';
options.ncpus = 4;
options.scope = 'null';
%%
[v2,adj2,e2] = gibbs(vdata,adj,edata, options);
%% display output

outputimg = zeros(imgdim);
for i = 1:imgdim
    for j = 1:imgdim
        outputimg(i,j) = v2{(i-1)*imgdim+j}.sample;
    end
end
figure;
imagesc(outputimg);
colormap('gray');
