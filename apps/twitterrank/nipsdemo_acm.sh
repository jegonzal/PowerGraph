#!/bin/sh
GLHOME="$HOME/graphlab2.2"
GRAPHHOME="$HOME/data/acmgraph"
pagerank_only="false"  #only run pagerank
lda_only="true"  #only run lda

cmd="$GLHOME/release/apps/twitterrank/twitterrank --lda_edges $GRAPHHOME/articlewords --lda_dictionary $GRAPHHOME/words --join_on_id=true --wordid_offset=1 --ncpus=8 --ntopics=20 --topklda=30 --pagerank_only $pagerank_only --lda_only $lda_only"

mpicmd="mpiexec -np 12 -hostfile $HOME/machines $cmd"
echo $mpicmd
exec $mpicmd
