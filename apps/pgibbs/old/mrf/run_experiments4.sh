#!/bin/bash

homedir="/mnt/bigbrofs/usr3/jegonzal"
root="$homedir/graphlab/version1/release/apps/pgibbs"
binary="pgibbs_denoise"
problem="$homedir/denoise_experiments/checkerboard3/problem.bin"


time_vals="{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,25,30,35,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000,1500,2000}"


printf "\n\n\n"
command="$root/$binary
            --problem=$problem 
            --experiment=colored 
            --ncpus=16
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
            --ncpus=16
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
            --ncpus=16
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
            --ncpus=16
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
            --ncpus=16
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
            --ncpus=16
            --times=$time_vals
            --weights=true
            --pruning=1E-1
            --nroots=16"
echo $command
printf "\n"
date
printf "\n"
time eval $command




