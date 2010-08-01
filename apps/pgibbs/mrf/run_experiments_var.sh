#!/bin/bash

homedir="/mnt/bigbrofs/usr3/jegonzal"
root="$homedir/graphlab/version1/release/apps/pgibbs"
binary="pgibbs_denoise"
problem="$homedir/denoise_experiments/checkerboard3/problem.bin"


time_vals="{2,4,8,10,15,20,25}"

for i in `seq 1 30`;
do

printf "\n\n\n"
command="$root/$binary
            --problem=$problem 
            --experiment=colored 
            --ncpus=1
            --times=$time_vals"
echo $command
printf "\n"
date
printf "\n"
time eval $command


printf "\n\n\n"
command="$root/$binary
            --problem=$problem 
            --experiment=async 
            --scheduler=multiqueue_fifo 
            --ncpus=1
            --times=$time_vals"
echo $command
printf "\n"
date
printf "\n"
time eval $command


printf "\n\n\n"
command="$root/$binary
            --problem=$problem
            --experiment=tree
            --scheduler=multiqueue_fifo
            --ncpus=1
            --times=$time_vals
            --nroots=16"
echo $command
printf "\n"
date
printf "\n"
time eval $command


printf "\n\n\n"
command="$root/$binary
            --problem=$problem
            --experiment=tree
            --treeheight=50
            --scheduler=multiqueue_fifo
            --ncpus=1
            --times=$time_vals
            --nroots=16"
echo $command
printf "\n"
date
printf "\n"
time eval $command



printf "\n\n\n"
command="$root/$binary
            --problem=$problem
            --experiment=tree
            --scheduler=multiqueue_priority
            --ncpus=1
            --times=$time_vals
            --weights=true
            --nroots=16"
echo $command
printf "\n"
date
printf "\n"
time eval $command


printf "\n\n\n"
command="$root/$binary
            --problem=$problem
            --experiment=tree
            --scheduler=multiqueue_priority
            --ncpus=1
            --times=$time_vals
            --weights=true
            --pruning=1E-1
            --nroots=16"
echo $command
printf "\n"
date
printf "\n"
time eval $command



done
