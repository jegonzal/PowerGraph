#!/bin/sh
GLHOME="/home/haijieg/graphlabapi"
GRAPHHOME="/home/haijieg/acmgraph"
pagerank_only="false"  #only run pagerank
lda_only="false"  #only run lda
personal_weights="./weights" #input file of the personalized topic weights

cmd="$GLHOME/debug/apps/twitterrank/twitterrank --pagerank_edges $GRAPHHOME/articlecitation --lda_edges $GRAPHHOME/articlewords --lda_dictionary $GRAPHHOME/words --join_on_id=true --wordid_offset=1 --ncpus=4 --ntopics=10 --topkpr=10 --topklda=20 --default_ndocs=1 --default_topicval=1 --pagerank_only $pagerank_only --lda_only $lda_only"
#--w_personal $personal_weights"

mpicmd="mpiexec -np 8 -hostfile /home/haijieg/machines $cmd"
echo $mpicmd
exec $mpicmd
