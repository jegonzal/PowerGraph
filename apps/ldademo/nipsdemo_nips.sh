#!/bin/sh
GLHOME="$HOME/graphlabapi"
GRAPHHOME="$HOME/data/nips"

cmd="$GLHOME/release/apps/ldademo/ldademo --lda_edges $GRAPHHOME/nipstext.txt --lda_dictionary $GRAPHHOME/nipswordlist.txt --wordid_offset=0 --ntopics=20 --topklda=30"

mpicmd="mpiexec -np 8 -hostfile $HOME/machines $cmd"
echo $mpicmd
exec $mpicmd
