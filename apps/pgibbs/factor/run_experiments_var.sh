#!/bin/bash

homedir="/mnt/bigbrofs/usr3/jegonzal"
root="$homedir/graphlab/version1/release/apps/pgibbs/factor"
binary="pgibbs_factor_alchemy"
network="$homedir/pgibbs_experiments/alchemy/cora-1/cora-1.alchemy"


time_vals="{1,2,4,8,10,15,20,25}"

for i in `seq 1 50`;
do

printf "\n\n\n"
command="$root/$binary
            --network=$network 
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
            --network=$network 
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
            --network=$network
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
            --network=$network
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
            --network=$network
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
            --network=$network
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



done
