clear;

D = 20;
nusers = 100;
nmovies = 1000;


ufactors = 2*randn(nusers, D);
mfactors = 2*randn(nmovies, D);

fid = fopen('synth_data.tsv','w');

for u = 1:nusers 
  movies = randperm(nmovies);
  for m = 1:100
    uid = u;
    mid = movies(m);
    rating = ufactors(uid) * mfactors(mid)' + randn()/2;
    fprintf(fid, '%d\t%d\t%f\n', uid, mid, rating);
  end
end

fclose(fid);
    