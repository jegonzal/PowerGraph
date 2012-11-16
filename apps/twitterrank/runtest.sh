#!/bin/sh
GLHOME="/home/haijieg/graphlabapi"
GRAPHHOME="/home/haijieg/test_graph"
cmd="$GLHOME/debug/apps/twitterrank/twitterrank --pagerank_edges $GRAPHHOME/pagerank/edges --pagerank_vertices $GRAPHHOME/pagerank/vertices --lda_edges $GRAPHHOME/lda/edges --lda_vertices $GRAPHHOME/lda/vertices --lda_dictionary $GRAPHHOME/lda/dictionary/part0 --join_on_id=false --wordid_offset=0 --ncpus=1 --ntopics=5 --topkpr=3 --topklda=5 --default_ndocs=1 --default_topicval=1"

# echo $cmd
# exec $cmd

mpicmd="mpiexec -np 2 -hostfile /home/haijieg/machines $cmd"
echo $mpicmd
exec $mpicmd
