clear;
raw = importdata('directed_triangles');
%%
col.vid       = 1;
col.in        = 3; % Triangles among people following you
col.out       = 2; % Triangles among peopel you follow
col.through   = 4;
col.cycle     = 5;
col.followers = 6;
col.following = 7;


%% compute top by cycle

degree = raw(:, col.followers) + raw(:, col.following);
[~,ind.degree] = sort(degree, 'descend');

[~,ind.followers] = sort(raw(:,col.followers), 'descend');

[~,ind.cycle] = sort(raw(:,col.cycle), 'descend');
[~,ind.in] = sort(raw(:,col.in), 'descend');
[~,ind.out] = sort(raw(:,col.out), 'descend');
[~,ind.through] = sort(raw(:,col.through), 'descend');

%%
cluster_coeff = sum(raw(:, [2,3,4,5]),2) ./ (degree + 1);
[~,ind.cluster_coeff] = sort(cluster_coeff, 'descend');

%%

cluster_coeff2 = sum(raw(:, [2,3,4,5]),2) ./ (raw(:, col.following) +1);
[~,ind.cluster_coeff2] = sort(cluster_coeff2, 'descend');

%% 

cluster_coeff3 = raw(:, col.through) ./ (degree + 1);
[~,ind.cluster_coeff3] = sort(cluster_coeff3, 'descend');



%% 
cluster_coeff4 = raw(:, col.in) ./ (raw(:, col.followers) + 1);
[~,ind.cluster_coeff4] = sort(cluster_coeff4, 'descend');


%% Render json
nusers = 10;


fid = fopen('top_users.json', 'w');
fprintf(fid, '[\n');

fprintf(fid, '\t { "name": "degree", "label": "Degree", "values": [\n');
for i = 1:nusers
   sep = ',\n';
   if(i == nusers)
        sep = '\n';
   end
   fprintf(fid, ['\t\t ["%d", "%d"]', sep], ...
        [raw(ind.degree(i), col.vid), degree(ind.degree(i))]); 
end
fprintf(fid, '\t]},\n');

fprintf(fid, '\t { "name": "followers", "label": "Followers", "values": [\n');
for i = 1:nusers
   sep = ',\n';
   if(i == nusers)
        sep = '\n';
   end
   fprintf(fid, ['\t\t ["%d", "%d"]', sep], ...
        raw(ind.followers(i), [col.vid, col.followers])); 
end
fprintf(fid, '\t]},\n');


fprintf(fid, '\t { "name": "cycle", "label": "Cycle Triangles", "values": [\n');
for i = 1:nusers
   sep = ',\n';
   if(i == nusers)
        sep = '\n';
   end
   fprintf(fid, ['\t\t ["%d", "%d"]', sep], ...
        raw(ind.cycle(i), [col.vid, col.cycle])); 
end
fprintf(fid, '\t]},\n');

% fprintf(fid, '\t { "name": "in", "label": "In Triangles", "values": [\n');
% for i = 1:nusers
%    sep = ',\n';
%    if(i == nusers)
%         sep = '\n';
%    end
%    fprintf(fid, ['\t\t ["%d", "%d"]', sep], ...
%         raw(ind.in(i), [col.vid, col.in])); 
% end
% fprintf(fid, '\t]},\n');
% 
% 
% fprintf(fid, '\t { "name": "out", "label": "Out Triangles", "values": [\n');
% for i = 1:nusers
%    sep = ',\n';
%    if(i == nusers)
%         sep = '\n';
%    end
%    fprintf(fid, ['\t\t ["%d", "%d"]', sep], ...
%         raw(ind.out(i), [col.vid, col.out])); 
% end
% fprintf(fid, '\t]},\n');
% 
% 
% fprintf(fid, '\t { "name": "through", "label": "Through Triangles", "values": [\n');
% for i = 1:nusers
%    sep = ',\n';
%    if(i == nusers)
%         sep = '\n';
%    end
%    fprintf(fid, ['\t\t ["%d", "%d"]', sep], ...
%         raw(ind.through(i), [col.vid, col.through])); 
% end
% fprintf(fid, '\t]},\n');
% 

fprintf(fid, '\t { "name": "cluster", "label": "Triangles / Degree", "values": [\n');
for i = 1:nusers
   sep = ',\n';
   if(i == nusers)
        sep = '\n';
   end
   fprintf(fid, ['\t\t ["%d", "%d"]', sep], ...
        [raw(ind.cluster_coeff(i), col.vid), cluster_coeff(ind.cluster_coeff(i))]); 
end
fprintf(fid, '\t]},\n');


fprintf(fid, '\t { "name": "cluster2", "label": "Triangles / Following", "values": [\n');
for i = 1:nusers
   sep = ',\n';
   if(i == nusers)
        sep = '\n';
   end
   fprintf(fid, ['\t\t ["%d", "%d"]', sep], ...
        [raw(ind.cluster_coeff2(i), col.vid), cluster_coeff2(ind.cluster_coeff2(i))]); 
end
fprintf(fid, '\t]},\n');

fprintf(fid, '\t { "name": "cluster3", "label": "Through Triangles / Degree", "values": [\n');
for i = 1:nusers
   sep = ',\n';
   if(i == nusers)
        sep = '\n';
   end
   fprintf(fid, ['\t\t ["%d", "%d"]', sep], ...
        [raw(ind.cluster_coeff3(i), col.vid), cluster_coeff3(ind.cluster_coeff3(i))]); 
end
fprintf(fid, '\t]},\n');



fprintf(fid, '\t { "name": "cluster4", "label": "In Triangles / Followers", "values": [\n');
for i = 1:nusers
   sep = ',\n';
   if(i == nusers)
        sep = '\n';
   end
   fprintf(fid, ['\t\t ["%d", "%d"]', sep], ...
        [raw(ind.cluster_coeff4(i), col.vid), cluster_coeff4(ind.cluster_coeff4(i))]); 
end
fprintf(fid, '\t]}\n');



fprintf(fid, ']\n');
fclose(fid);
