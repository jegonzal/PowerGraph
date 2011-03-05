function asgs = spectral_bisection(E,levels)
%% Force E to be symmetric
M = E + E';

asgs = ones(length(E),1);
new_asgs = asgs;

for i = 1:levels
   for j = unique(asgs)'
      next_class = max(new_asgs) + 1;
      ind = asgs == j;
      Msub = M(ind, ind);
      degree = sum(Msub);
      L = diag(degree) - Msub;
      [V,D] = eigs(L, 2, 'sm');
      new_ind = V(:,1) > 0;
      revmap = find(ind);
      new_asgs(revmap(new_ind)) = next_class;
   end
   asgs = new_asgs;
end



end





