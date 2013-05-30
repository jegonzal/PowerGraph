home="/Users/haijieg"
glhome="/Users/haijieg/graphlab/graphlab2.2"
graph="$home/data/netflix/small/feat_netflix_edge.csv"
mvname="$home/data/netflix/meta/movienames.short.txt"
genre="$home/data/netflix/small/feat_netflix_genre.csv"
topics="$home/data/netflix/small/feat_netflix_topic.csv"
bin="$glhome/release/apps/netflix++/netflix_main"

# --genre_feature=$genre \
# --topic_feature=$topics \
mpiexec -np 2 $bin --matrix=$graph --D=20 \
  --movielist=$mvname \
  --use_bias=false \
  --use_als=true \
  --use_feature_weights=false \
  --use_feature_latent=false \
  --max_iter=5 \
  --interactive=true \
  --testpercent=0.2
