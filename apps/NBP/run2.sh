#!/bin/sh
#./NBP --kbpfile=/mnt/bigbrofs/usr5/lesong/data/make3ddata/KBP_fold100 --gmmfile=/mnt/bigbrofs/usr5/lesong/data/make3ddata/GMM_fold100 --iterations=2 --ncpus=6
DATAFILES=/mnt/bigbrofs/usr5/lesong/data/make3ddata/
uname -a
date
for i in {1..274}
do
./NBP --gmmfile=$DATAFILES/GMM_fold${i} --kbpfile=$DATAFILES/KBP_fold${i} --iterations=10 --ncpus=8  --epsilon=.1 | tee /mnt/bigbrofs/usr6/bickson/`date +%F%T.nbp3D_big`.IT${i}
date
done
