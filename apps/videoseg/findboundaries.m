function boundaries = findboundaries(seg)


tmp = zeros([3,3,3]);
tmp(2,2,2) = 6;

tmp(2,3,2) = -1;
tmp(2,1,2) = -1;
tmp(1,2,2) = -1;
tmp(3,2,2) = -1;
tmp(2,2,1) = -1;
tmp(2,2,3) = -1;

pos = imfilter(seg, tmp,'replicate');


boundaries = find(pos);

end