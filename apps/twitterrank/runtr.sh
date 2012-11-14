#!/bin/sh
GLHOME="/home/haijieg/graphlabapi"
GRAPHHOME="/home/haijieg/test_graph"
cmd="$GLHOME/release/apps/twitterrank/twitterrank --pagerank_edges $GRAPHHOME/pagerank/edges --pagerank_vertices $GRAPHHOME/pagerank/vertices --lda_edges $GRAPHHOME/lda/edges --lda_vertices $GRAPHHOME/lda/vertices --lda_dictionary $GRAPHHOME/lda/dictionary/part0 --ntopics 5"
echo $cmd
exec $cmd
