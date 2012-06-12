clear;

D = 20;
nusers = 10000;
nmovies = 100000;


ufactors = 2*randn(nusers, D);
mfactors = 2*randn(nmovies, D);


!mkdir synth_d20

fid_train = fopen('synth_d20/train_data.tsv','w');
fid_validate = fopen('synth_d20/train_data.tsv.validate','w');
fid_test = fopen('synth_d20/train_data.tsv.test','w');
for u = 1:nusers 
  movies = randperm(nmovies);
  for m = 1:50
    uid = u;
    mid = movies(m);
    rating = ufactors(uid) * mfactors(mid)' + randn()/2;
    fprintf(fid, '%d\t%d\t%f\n', uid, mid, rating);
  end
end
fclose(fid);
    