function nbrs = neighborsof(segmat, idx)

pos= zeros(size(segmat));
pos(segmat == idx) = 1;

tmp = zeros([3,3,3]);
tmp(2,2,2) = 1;
tmp(2,3,2) = 1;
tmp(2,1,2) = 1;
tmp(1,2,2) = 1;
tmp(3,2,2) = 1;
tmp(2,2,1) = 1;
tmp(2,2,3) = 1;
%strel(tmp);

pos = imdilate(pos, tmp);


nbrs = unique(segmat(find(pos)));
nbrs = nbrs(nbrs ~= idx);
end