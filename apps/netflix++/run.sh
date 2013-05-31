home="/Users/haijieg"
glhome="/Users/haijieg/graphlab/graphlab2.2"
graph="$home/data/netflix/small/feat_netflix_edge.csv"
mvname="$home/data/netflix/meta/moviename_escaped.txt"
genre="$home/data/netflix/small/feat_netflix_genre.csv"
topics="$home/data/netflix/small/feat_netflix_topic.csv"
bin="$glhome/release/apps/netflix++/netflix_main"

# --genre_feature=$genre \
# --topic_feature=$topics \
mpiexec -np 4 $bin --matrix=$graph --D=10 \
  --movielist=$mvname \
  --use_bias=false \
  --use_als=true \
  --use_feature_weights=false \
  --use_feature_latent=false \
  --max_iter=20 \
  --interactive=true \
  --testpercent=0.2
