#!/bin/sh
GLHOME="$HOME/graphlab2.2"
GRAPHHOME="$HOME/data/nips"
pagerank_only="false"  #only run pagerank
lda_only="true"  #only run lda

cmd="$GLHOME/release/apps/twitterrank/twitterrank --execution synchronous --lda_edges $GRAPHHOME/nipstext.txt --lda_dictionary $GRAPHHOME/nipswordlist.txt --wordid_offset=0 --ntopics=20 --topklda=30 --pagerank_only $pagerank_only --lda_only $lda_only"

mpicmd="mpiexec -np 12 -hostfile $HOME/machines $cmd"
# mpicmd="$cmd"
echo $mpicmd
exec $mpicmd
