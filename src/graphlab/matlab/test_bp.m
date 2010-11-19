arity = 5;
imgdim = 100;
%% generate images
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

% make noisy image
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
        vtx.belief = ones(1,arity);
        vtx.unary = pdf('norm',0:(arity-1), noisyimg(i,j), 1);
        vtx.unary = vtx.unary ./ sum(vtx.unary);
        vdata{(i-1)*100+j} = vtx;
    end
end

%% generate ising edgedata
edge.binary = zeros(arity);
for i = 1:arity
    for j = 1:arity
        edge.binary(i,j) = exp(-3 * abs(i-j));
    end
end
edge.msg = ones(1,5);
edata = {edge};
%% generate adjacency matrix
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
% you need to set the graphlab directory
compile_update_function({'bp_update'}, vdata{1},edata{1}, [getenv('HOME') '/graphlabapi/graphlabapi/release/src/graphlab'], 'bp', 'bp');
%% set options
options.initial_schedule(1).update_function = 'bp_update';
options.initial_schedule(1).vertices=uint32(1:(imgdim * imgdim));
options.initial_schedule(1).priorities=ones(size(options.sched.vertices));
options.scheduler = 'sweep';
options.ncpus = 4;
options.scope = 'edge';
%% cd in the BP directory and call
cd bp
[v2,adj2,e2] = bp(vdata,adj,edata, options);
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
