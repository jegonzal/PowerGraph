function [chunks, groupname, cuts] = group(M, group_on, aggregate)
% Example Usage:
% grps = group(raw.data, [1,2], [1,2,4]);
% groups together (puts in the same cell) the rows of raw.data that have
% the same values for columns [1,2]. And the data that is put in the cell
% is the values in cols [1,2,4].


%% This is a very powerful function to group large tables
%  Let M be a large table with each row being a record
%  then group on are the columns on which to group the records
%  and aggregate are the columns to retain in the final grouping
%  chunks will be a cell array where each cell contains all the
%  records with groupname and cuts is the count of the items in that
%  group.  
[groupname, junk, ind] = unique(M(:,group_on), 'rows');
clear('junk');
[junk, oind] = sort(ind);
clear('junk');
cuts = hist(ind, min(ind):max(ind))';
chunks = mat2cell(M(oind,aggregate), cuts, length(aggregate));
end