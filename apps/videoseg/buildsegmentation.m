function segmatrix = buildsegmentation(vidmat,numpart)
numpart = numpart - 1;
% build a graph
numframes = size(vidmat, 1);
height = size(vidmat, 2);
width = size(vidmat, 3);
segmatrix = zeros([numframes,height,width]);
%segmatrix = vidmat(:,:,:,1)/2;
% compute all the pairwise pixel differences
Xdiff = zeros(numframes, height, width);
Ydiff = zeros(numframes, height, width);

parfor i = 1:numframes
Xdiff(i,:,:) = round(abs(imfilter(squeeze(vidmat(i,:,:)),fspecial('sobel'))) ) ;
Ydiff(i,:,:) = round(abs(imfilter(squeeze(vidmat(i,:,:)),fspecial('sobel')')) ) ;
end
d = max([max(max(max(Xdiff))),max(max(max(Ydiff)))]);
%Xdiff = (d + 1) - Xdiff;
%Ydiff = (d + 1) - Ydiff;
Xdiff = 1+100*exp(-Xdiff / d);
Ydiff = 1+100*exp(-Ydiff / d);
interframeC = 2;
opposingC = 25;

display('built edge images');
[adj, xyz] = grid_graph(numframes, height, width);
% this includes self edges. Remove them
% build the graph
adj = adj - speye(size(adj));
[row, col] = find(adj);

clear adj
pix1coor = xyz(row,:);
pix2coor = xyz(col,:);


valx=ones(size(row));

for i = 1:length(row)
    p = min(pix1coor(i,:), pix2coor(i,:));
    if (pix1coor(i,1) ~= pix2coor(i,1))
        valx(i)=interframeC;
    elseif (pix1coor(i,2) ~= pix2coor(i,2))
        valx(i) = Xdiff(p(1),p(2),p(3));
    elseif (pix1coor(i,3) ~= pix2coor(i,3))
        valx(i) = opposingC;
    end
    if (mod(i, 1000000) == 0)
        disp([num2str(i), '/', num2str(length(row))]);
    end
end
% clear some stuff I no longer need
clear Xdiff; 
%adj = sparse(row,col,val);

display('constructed adjacency matrix');

% add vertex s and t
% cut across frames
infinity = 10000000;
s = numframes*height*width+1;
t = numframes*height*width+2;



% cut across height
firstframeidx = find(xyz(:,2) == 1);
lastframeidx = find(xyz(:,2) == height);
numffirst = length(firstframeidx);
numflast = length(lastframeidx);

adj2 = sparse([row; s*ones(numffirst,1) ; lastframeidx],  ...
              [col;firstframeidx ; t*ones(numflast,1)],   ...
               [valx; infinity*ones(numffirst+numflast,1)]);
clear valx
adj2(t,t) = 0;
display 'cutting height'
rpart = rand(1,numpart);
for i = 1:numpart
    display(num2str(i))
    [maxflow, mincut,R,F] = max_flow(adj2, s, t);
    maxflow
    mincut = mincut(1:(s-1));
    edges = (mincut(row) == 1 & mincut(col) == -1) | (mincut(row) == -1 & mincut(col) == 1);
    % set the cut to have infinite weight
    % row(edges) are vertices on the boundary
    % I want every boundary vertex <=> boundary vertex edge to have weight
    % infinity
    u = unique(row(edges));
    u = dilateelem(size(segmatrix),u);
    selector = ismember(row,u) & ismember(col,u);
    adj3=sparse(row(selector),col(selector),infinity);
    adj3(t,t)=0;
    adj2 = adj2 + adj3;
    
    %adj3=sparse(row(edges),col(edges),infinity);
    %adj3(t,t)=0;
    %adj2 = adj2 + adj3;
    segmatrix(mincut == 1) = segmatrix(mincut == 1) + rpart(i);
end

[u,~,n] = unique(segmatrix(:));
u=1:length(u);
segmatrix=reshape(u(n),[numframes,height,width]);
length(u)








valy=ones(size(row));

for i = 1:length(row)
    p = min(pix1coor(i,:), pix2coor(i,:));
    if (pix1coor(i,1) ~= pix2coor(i,1))
        valy(i)=interframeC;
    elseif (pix1coor(i,2) ~= pix2coor(i,2))
        valy(i) = opposingC;
    elseif (pix1coor(i,3) ~= pix2coor(i,3))
        valy(i) = Ydiff(p(1),p(2),p(3));
    end
    if (mod(i, 1000000) == 0)
        disp([num2str(i), '/', num2str(length(row))]);
    end
end

clear Ydiff;

% cut across width
firstframeidx = find(xyz(:,3) == 1);
lastframeidx = find(xyz(:,3) == width);
numffirst = length(firstframeidx);
numflast = length(lastframeidx);

adj2 = sparse([row; s*ones(numffirst,1) ; lastframeidx],  ...
              [col;firstframeidx ; t*ones(numflast,1)],   ...
               [valy; infinity*ones(numffirst+numflast,1)]);
clear valy;
adj2(t,t) = 0;
display 'cutting width'
rpart = rand(1,numpart);
for i = 1:numpart
    display(num2str(i))
    [maxflow, mincut] = max_flow(adj2, s, t);
    maxflow
    mincut = mincut(1:(s-1));
    edges = (mincut(row) == 1 & mincut(col) == -1) | (mincut(row) == -1 & mincut(col) == 1);
    % set the cut to have infinite weight
    % row(edges) are vertices on the boundary
    % I want every boundary vertex <=> boundary vertex edge to have weight
    % infinity
    u = unique(row(edges));
    u = dilateelem(size(segmatrix),u);
    selector = ismember(row,u) & ismember(col,u);
    adj3=sparse(row(selector),col(selector),infinity);
    adj3(t,t)=0;
    adj2 = adj2 + adj3;
    %adj3=sparse(row(edges),col(edges),infinity);
    %adj3(t,t)=0;
    %adj2 = adj2 + adj3;
    segmatrix(mincut == 1) = segmatrix(mincut == 1) + rpart(i);
end

[u,~,n] = unique(segmatrix(:));
u=1:length(u);
segmatrix=reshape(u(n),[numframes,height,width]);
length(u)

end