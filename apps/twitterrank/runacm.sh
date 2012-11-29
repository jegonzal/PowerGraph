#!/bin/sh
GLHOME="$HOME/graphlabapi"
GRAPHHOME="$HOME/acmgraph"
pagerank_only="false"  #only run pagerank
lda_only="false"  #only run lda
personal_weights="./weights" #input file of the personalized topic weights

cmd="$GLHOME/release/apps/twitterrank/twitterrank --pagerank_edges $GRAPHHOME/articlecitation --lda_edges $GRAPHHOME/articlewords --lda_dictionary $GRAPHHOME/words --join_on_id=true --wordid_offset=1 --ncpus=8 --ntopics=20 --topkpr=10 --topklda=30 --default_ndocs=1 --default_topicval=1 --pagerank_only $pagerank_only --lda_only $lda_only --w_personal $personal_weights"

mpicmd="mpiexec -np 8 -hostfile $HOME/machines $cmd"
#mpicmd="$cmd"
echo $mpicmd
exec $mpicmd
