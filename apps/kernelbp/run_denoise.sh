#!/bin/bash
datadir="/mnt/bigbrofs/usr5/lesong/codes/make3ddata/images"

for i in 10 25 50 75 100 125 150 175 200 225 250
do

for t in {1..50}
do
./kbp ${datadir}/KBP_model${i}.it ${datadir}/Test${i}_${t}.it log.txt 8
./discrete_bp_itgmm --gmmfile=${datadir}/GMM_model${i}.it --inputfile=${datadir}/Test${i}_${t}.it --iterations=30 --ncpus=8 --scheduler=sweep --logfile=log.txt
./particle_bp_itgmm --gmmfile=${datadir}/GMM_model${i}.it --inputfile=${datadir}/Test${i}_${t}.it --particles=50 --iterations=30 --ncpus=8 --scheduler=sweep --resample=100 --logfile=log.txt
done

done