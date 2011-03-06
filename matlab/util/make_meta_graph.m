function [meta_eW, meta_vW] = make_meta_graph(eW, vW, cut)
[n,n] = size(vW);
[u,v,w] = find(eW);

cutw = w .* (cut(u) ~= cut(v));
cutu = cut(u);
cutv = cut(v);
parts = unique(cut(:));
meta_eW = zeros(max(cut));

[eweights, edges] = group([cutu, cutv, cutw], [1,2],[3]);
eweights = cellfun(@sum, eweights);
meta_eW = sparse(edges(:,1), edges(:,2), eweights, max(parts), max(parts));
meta_eW = full(meta_eW);

% for ui = parts';
%   disp(['ui: ', num2str(ui)]);
%   for uj = parts';
%     disp(['uj: ', num2str(uj)]);  
%     ind = (cutu == ui & cutv == uj);
%     meta_eW(ui,uj) = sum(cutw(ind));
%   end
% end

[vweights, parts] = group([cut, vW], 1, 2);
meta_vW = zeros(max(parts), 1);
meta_vW(parts) = cellfun(@sum, vweights);






end