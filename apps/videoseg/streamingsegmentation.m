function segmatrix = streamingsegmentation(vidmat,numpart, windowsize)
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
Xdiff(i,:,:) = round(abs(imfilter(squeeze(vidmat(i,:,:)),fspecial('sobel'),'circular')) ) ;
Ydiff(i,:,:) = round(abs(imfilter(squeeze(vidmat(i,:,:)),fspecial('sobel')','circular')) ) ;
end
d = max([max(max(max(Xdiff))),max(max(max(Ydiff)))]);
%Xdiff = (d + 1) - Xdiff;
%Ydiff = (d + 1) - Ydiff;
Xdiff = 1+100*exp(-Xdiff / d);
Ydiff = 1+100*exp(-Ydiff / d);
interframeC = 5;
opposingC = 25;

display('built edge images');
[adj, xyz] = grid_graph(windowsize, height, width);
% this includes self edges. Remove them
% build the graph
adj = adj - speye(size(adj));
[row, col] = find(adj);
clear adj
pix1coor = xyz(row,:);
pix2coor = xyz(col,:);


display('constructed adjacency matrix');

% add vertex s and t
% cut across frames
infinity = 10000000;
s = windowsize*height*width+1;
t = windowsize*height*width+2;
valx=ones(size(row));


for partnum = 1:numpart
    display(['partitioning: ',num2str(partnum)]);
    r = rand();
    lastframecuts = [];
    prevsegmatrix = segmatrix;
    for step = 1:windowsize-1:(numframes-windowsize+1)
        curframe = step;
        lastframe = step+windowsize-1;
        adds = [];
        addt = [];
        if (step == 1)
            for i = 1:length(row)
                p = min(pix1coor(i,:), pix2coor(i,:));
                if (pix1coor(i,1) ~= pix2coor(i,1))
                    valx(i)=interframeC;
                elseif (pix1coor(i,2) ~= pix2coor(i,2))
                    valx(i) = Xdiff(curframe + p(1) - 1,p(2),p(3));
                elseif (pix1coor(i,3) ~= pix2coor(i,3))
                    valx(i) = opposingC;
                end
                if (mod(i, 1000000) == 0)
                    disp([num2str(i), '/', num2str(length(row))]);
                end
            end
        else        
            % connect the cuts we made in the last frame
            for i = 1:length(row)
                p = min(pix1coor(i,:), pix2coor(i,:));
                % first frame
                if (pix1coor(i,1) == 1 && pix2coor(i,1) == 1)
                    x1=pix1coor(i,2); y1=pix1coor(i,3);
                    x2=pix2coor(i,2); y2=pix2coor(i,3);
                    
                    if (lastframecuts(x1,y1) ~= lastframecuts(x2,y2))
                        % boundary    
                        valx(i) = 0;
                        % remember this. we need to connect st to it
                        if (lastframecuts(x1,y1) == 1)
                            adds = [adds,row(i)];
                        elseif (lastframecuts(x1,y1) == -1)
                            addt = [addt,row(i)];
                        end
                        continue;
                    end
                end
                
                if (pix1coor(i,1) ~= pix2coor(i,1))
                    valx(i)=interframeC;
                elseif (pix1coor(i,2) ~= pix2coor(i,2))
                    valx(i) = Xdiff(curframe + p(1) - 1,p(2),p(3));
                elseif (pix1coor(i,3) ~= pix2coor(i,3))
                    valx(i) = opposingC;
                end
                if (mod(i, 1000000) == 0)
                    disp([num2str(i), '/', num2str(length(row))]);
                end
            end

        end
        adds=unique(adds)';
        addt=unique(addt)';
        if (numel(intersect(adds,addt)) > 0)
            display ('s matching t!!');
            pause
        end
        firstframeidx = find(xyz(:,2) == 1);
        lastframeidx = find(xyz(:,2) == height);
        numffirst = length(firstframeidx);
        numflast = length(lastframeidx);

        adj2 = sparse([row; s*ones(numffirst+length(adds),1) ; lastframeidx;addt],  ...
                      [col;firstframeidx;adds; t*ones(numflast+length(addt),1)],   ...
                       [valx; infinity*ones(numffirst+numflast+length(adds)+length(addt),1)]);
        adj2(t,t) = 0;
        
        % fill in the other cuts we have already made
        if (partnum > 1)
            curpart = prevsegmatrix(curframe:lastframe,:,:);
            u = findboundaries(curpart);
            u = dilateelem(size(curpart),u);
            selector = ismember(row,u) & ismember(col,u);
            selector = selector & (valx > 0); % do not add to edges which we set to 0
            adj3=sparse(row(selector),col(selector),infinity);
            adj3(t,t)=0;
            adj2 = adj2 + adj3;
        end


        [maxflow, mincut,R,F] = max_flow(adj2, s, t);
        mincut = mincut(1:(s-1));
        lastframecuts = reshape(mincut,[windowsize,height,width]);
        lastframecuts = squeeze(lastframecuts(windowsize,:,:));
        
        maxflow

        %adj3=sparse(row(edges),col(edges),infinity);
        %adj3(t,t)=0;
        %adj2 = adj2 + adj3;
        tmpmatrix = zeros([windowsize,height,width]);
        tmpmatrix(mincut == 1) = tmpmatrix(mincut == 1) + r;
        if (curframe > 1)
            segmatrix((curframe+1):lastframe,:,:) = segmatrix((curframe+1):lastframe,:,:) + tmpmatrix(2:end,:,:);
        else
            segmatrix((curframe):lastframe,:,:) = segmatrix((curframe):lastframe,:,:) + tmpmatrix;
        end
    end
end

[u,~,n] = unique(segmatrix(:));
u=1:length(u);
xsegmatrix=reshape(u(n),[numframes,height,width]);
length(u)












clear valx;







segmatrix = zeros([numframes,height,width]);

valy=ones(size(row));

% other direction
for partnum = 1:numpart
    display(['partitioning: ',num2str(partnum)]);
    r = rand();
    lastframecuts = [];
    prevsegmatrix = segmatrix;
    for step = 1:windowsize-1:(numframes-windowsize+1)
        curframe = step;
        lastframe = step+windowsize-1;
        adds = [];
        addt = [];
        if (step == 1)
            for i = 1:length(row)
                p = min(pix1coor(i,:), pix2coor(i,:));
                if (pix1coor(i,1) ~= pix2coor(i,1))
                    valy(i)=interframeC;
                elseif (pix1coor(i,2) ~= pix2coor(i,2))
                    valy(i) = opposingC;
                elseif (pix1coor(i,3) ~= pix2coor(i,3))
                    valy(i) = Ydiff(curframe + p(1) - 1,p(2),p(3));
                end
                if (mod(i, 1000000) == 0)
                    disp([num2str(i), '/', num2str(length(row))]);
                end
            end
        else        
            % connect the cuts we made in the last frame
            for i = 1:length(row)
                p = min(pix1coor(i,:), pix2coor(i,:));
                % first frame
                if (pix1coor(i,1) == 1 && pix2coor(i,1) == 1)
                    x1=pix1coor(i,2); y1=pix1coor(i,3);
                    x2=pix2coor(i,2); y2=pix2coor(i,3);
                    
                    if (lastframecuts(x1,y1) ~= lastframecuts(x2,y2))
                        % boundary    
                        valy(i) = 0;
                        % remember this. we need to connect st to it
                        if (lastframecuts(x1,y1) == 1)
                            adds = [adds,row(i)];
                        elseif (lastframecuts(x1,y1) == -1)
                            addt = [addt,row(i)];
                        end
                        continue;
                    end
                end
                
                if (pix1coor(i,1) ~= pix2coor(i,1))
                    valy(i)=interframeC;
                elseif (pix1coor(i,2) ~= pix2coor(i,2))
                    valy(i) = opposingC;
                elseif (pix1coor(i,3) ~= pix2coor(i,3))
                    valy(i) = Ydiff(curframe + p(1) - 1,p(2),p(3));
                end
                if (mod(i, 1000000) == 0)
                    disp([num2str(i), '/', num2str(length(row))]);
                end
            end

        end
        adds=unique(adds)';
        addt=unique(addt)';
        if (numel(intersect(adds,addt)) > 0)
            display ('s matching t!!');
            pause
        end
        firstframeidx = find(xyz(:,3) == 1);
        lastframeidx = find(xyz(:,3) == width);
        numffirst = length(firstframeidx);
        numflast = length(lastframeidx);

        adj2 = sparse([row; s*ones(numffirst+length(adds),1) ; lastframeidx;addt],  ...
                      [col;firstframeidx;adds; t*ones(numflast+length(addt),1)],   ...
                       [valy; infinity*ones(numffirst+numflast+length(adds)+length(addt),1)]);
        adj2(t,t) = 0;
        
        % fill in the other cuts we have already made
        if (partnum > 1)
            curpart = prevsegmatrix(curframe:lastframe,:,:);
            u = findboundaries(curpart);
            u = dilateelem(size(curpart),u);
            selector = ismember(row,u) & ismember(col,u);
            selector = selector & (valy > 0); % do not add to edges which we set to 0
            adj3=sparse(row(selector),col(selector),infinity);
            adj3(t,t)=0;
            adj2 = adj2 + adj3;
        end


        [maxflow, mincut,R,F] = max_flow(adj2, s, t);
        mincut = mincut(1:(s-1));
        lastframecuts = reshape(mincut,[windowsize,height,width]);
        lastframecuts = squeeze(lastframecuts(windowsize,:,:));
        
        maxflow

        %adj3=sparse(row(edges),col(edges),infinity);
        %adj3(t,t)=0;
        %adj2 = adj2 + adj3;
        tmpmatrix = zeros([windowsize,height,width]);
        tmpmatrix(mincut == 1) = tmpmatrix(mincut == 1) + r;
        if (curframe > 1)
            segmatrix((curframe+1):lastframe,:,:) = segmatrix((curframe+1):lastframe,:,:) + tmpmatrix(2:end,:,:);
        else
            segmatrix((curframe):lastframe,:,:) = segmatrix((curframe):lastframe,:,:) + tmpmatrix;
        end
    end
end

segmatrix = segmatrix + xsegmatrix;
[u,~,n] = unique(segmatrix(:));
u=1:length(u);
segmatrix=reshape(u(n),[numframes,height,width]);
length(u)






end