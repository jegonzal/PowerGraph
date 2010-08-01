function m = videotomat(vid, numframes)
numframes=min([numframes,vid.NumberOfFrames]);
m = zeros([numframes, vid.Height, vid.Width]);

for i = 1:numframes
    m(i,:,:) = imfilter(rgb2gray(read(vid,i)),fspecial('gaussian'));
end