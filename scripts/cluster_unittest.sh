#/bin/sh
cd /tmp/graphlabapi/release/demoapps/pmf
ln -s ~/www-graphl1ab/graphlab/datasets/movielens .
echo "Testing LDA+clusterdump"
./glcluster movielens 3 10 0 --float=true --pmfformat=true --max_iter=10 --clusterdump=true
echo "testing KMEANS"
./glcluster movielens 0 10 0 --float=true --pmfformat=true --max_iter=10 --binaryoutput=true
echo "testing Fuzzy Kmeans"
./glcluster movielens 2 10 0 --float=true --pmfformat=true --max_iter=10 --binaryoutput=true
