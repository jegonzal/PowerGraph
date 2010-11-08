#!/bin/bash 
datadir="/mnt/bigbrofs/usr5/lesong/codes/make3ddata/images"
uname -a
date
for t in {1..10}
do
for i in 10 25 50 75 100 125 150 175 200 225 250
do
  ./NBPD  --gmmfile=${datadir}/GMM_model${i}.it --inputfile=${datadir}/Test${i}_${t}.it  --iterations=30 --ncpus=8 --epsilon=1 |  tee /mnt/bigbrofs/usr6/bickson/`date +%F%T.nbp_denoise`.it${i}.${t}

date
done
done
