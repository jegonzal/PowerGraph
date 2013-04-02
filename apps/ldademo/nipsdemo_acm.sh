#!/bin/sh
GLHOME="$HOME/graphlabapi"
GRAPHHOME="$HOME/data/acmgraph"

cmd="$GLHOME/release/apps/ldademo/ldademo --lda_edges $GRAPHHOME/articlewords --lda_dictionary $GRAPHHOME/words --wordid_offset=1 --ntopics=20 --topklda=30"

mpicmd="mpiexec -np 8 -hostfile $HOME/machines $cmd"
echo $mpicmd
exec $mpicmd
