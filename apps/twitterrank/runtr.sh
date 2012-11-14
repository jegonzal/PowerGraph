#!/bin/sh
GLHOME="/home/haijieg/graphlabapi"
cmd="$GLHOME/release/apps/twitterrank/twitterrank --pagerank_edges /Users/haijieg/tweets/test_graph/pagerank/edges --pagerank_vertices /Users/haijieg/tweets/test_graph/pagerank/vertices --lda_edges /Users/haijieg/tweets/test_graph/lda/edges --lda_vertices /Users/haijieg/tweets/test_graph/lda/vertices --lda_dictionary /Users/haijieg/tweets/test_graph/lda/dictionary --ntopics 5"
echo $cmd
exec $cmd
