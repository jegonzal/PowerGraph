function score = score_cut(wE, wV, cut) 
[n,n] = size(wE);
[u,v,w] = find(wE);

p = length(unique(cut));

vCut = full(sum(...
   sparse(u,v, w .* (cut(u) ~= cut(v)), n,n) ...
   ))';
degW = sum(wE);
 
% For each processor compute the work and communication done by tha
% processor
nparts = max(unique(cut));
work = zeros(nparts,1);
ework = zeros(nparts,1);
comm = zeros(nparts,1);

parfor i = 1:nparts;
  %disp(['Scoring cut: ', num2str(i)]);
  ind = logical(cut == i);
  work(i) = sum( wV(ind) );
  ework(i) = sum( degW(ind) );
  comm(i) = sum( vCut(ind) );
end

score.work = work;
score.ework = ework;
score.comm = comm;
score.cost = sum(vCut(:)); %sum(score.comm);
score.bal = max(score.work) * p / sum(wV);


score.ebal = max(score.ework) * p / sum(wE(:));

% score.value = max(score.work) + sum(vCut);
end