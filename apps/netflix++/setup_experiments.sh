home="/Users/haijieg"
glhome="/Users/haijieg/graphlab/graphlab2.2"
graph="$home/data/netflix/small/splits"
mvname="$home/data/netflix/meta/moviename_escaped.txt"
genre="$home/data/netflix/small/feat_netflix_genre.csv"
topics="$home/data/netflix/small/feat_netflix_topic.csv"
bin="$glhome/release/apps/netflix++/netflix_main"
np=4

outdirbase="$glhome/apps/netflix++/output"

# nlatent_list=( 5 10 20 50 100 )
# iter_list=( 4 8 12 16 20 )
nlatent_list=( 5 )
iter_list=( 4 )
lambda_list=( 0.01 )

if [ -d "$outdirbase" ]
then
  mv "$outdirbase" "$outdirbase.orig"
fi

for NLATENT in "${nlatent_list[@]}"
do
  for ITER in "${iter_list[@]}"
  do
    for LAMBDA in "${lambda_list[@]}"
    do
      outdir="$outdirbase/D=${NLATENT}_iter=${ITER}_lambda=${LAMBDA}";
      mkdir -p "$outdir"
      cmd="mpiexec -np $np $bin --matrix=$graph \
        --D=$NLATENT \
        --lambda=$LAMBDA \
        --movielist=$mvname \
        --use_bias=false \
        --use_als=true \
        --use_feature_weights=false \
        --use_feature_latent=false \
        --max_iter=$ITER \
        --interactive=false \
        --saveprefix="$outdir/result"
        --testpercent=-1"
      echo $cmd > "$outdir/run.sh"
      chmod u+x "$outdir/run.sh"
    done
  done
done
echo 'for dir in D=*; do $dir/run.sh; done' > "$outdirbase/runall.sh"
chmod u+x "$outdirbase/runall.sh"
