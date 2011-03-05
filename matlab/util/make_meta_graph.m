function [meta_eW, meta_vW] = make_meta_graph(eW, vW, cut)
[n,n] = size(vW);
[u,v,w] = find(eW);

cutw = w .* (cut(u) ~= cut(v));
cutu = cut(u);
cutv = cut(v);
parts = unique(cut(:));
meta_eW = zeros(max(cut));

[weights, edges] = group([cutu, cutv, cutw], [1,2],[3]);
weights = cellfun(@sum, weights);
meta_eW = sparse(edges(:,1), edges(:,2), weights, max(parts), max(parts));
meta_eW = full(meta_eW);

% for ui = parts';
%   disp(['ui: ', num2str(ui)]);
%   for uj = parts';
%     disp(['uj: ', num2str(uj)]);  
%     ind = (cutu == ui & cutv == uj);
%     meta_eW(ui,uj) = sum(cutw(ind));
%   end
% end

meta_vW = 0;






end