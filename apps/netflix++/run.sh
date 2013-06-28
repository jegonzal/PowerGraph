home="/home/haijieg"
glhome="$home/graphlab/graphlab-dev"
graph="/data/netflix/features/joined_edge_list.csv"
mvname="/data/netflix/meta/moviename.txt"
genre="/data/netflix/features/feat_netflix_genre.csv"
topics="/data/netflix/features/feat_netflix_topic.csv"
bin="$glhome/release/apps/netflix++/netflix_main"

# --genre_feature=$genre \
# --topic_feature=$topics \

mpiexec -np 4 $bin --matrix=$graph --D=5 \
  --movielist=$mvname \
  --use_bias=true \
  --use_als=true \
  --use_feature_weights=true \
  --use_feature_latent=true \
  --max_iter=1 \
  --interactive=true \
  --testpercent=0.2 \
  --saveprefix="output/result"
