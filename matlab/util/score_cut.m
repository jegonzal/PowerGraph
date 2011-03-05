function score = score_cut(cut, wE, wV) 
[n,n] = size(wE);
[u,v,w] = find(wE);

p = length(unique(cut));

vCut = full(sum(...
   sparse(u,v, w .* (cut(u) ~= cut(v)), n,n) ...
   ))';

% For each processor compute the work and communication done by tha
% processor
for i = unique(cut(:))';
   score.work(i) = sum( wV(cut == i) );
   score.comm(i) = sum( vCut(cut == i) );
end

score.cost = sum(score.comm);
score.bal = max(score.work) * p / sum(wV);

% score.value = max(score.work) + sum(vCut);
end