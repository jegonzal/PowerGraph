#!/bin/bash

homedir="/mnt/bigbrofs/usr3/jegonzal"
root="$homedir/graphlab/version1/release/apps/pgibbs"
binary="pgibbs_denoise"
problem="$homedir/denoise_data/hard3/problem.bin"

#cpu_seq="1 `seq 2 2 16`"
cpu_seq=16

samples="1 2 3 4 5 6 7 8 9 10 15 20 30 40 50 100 150 250 350 450 500 550 600 650 700 1000 1500 2000 5000 10000 50000"


for nsamples in $samples; do
    for ncpus in $cpu_seq; do
        echo "###################################################################"
        echo "Using ncpus: $ncpus"
        echo "using samples: $nsamples"
        echo "-------------------------------------------------------------------"
        echo "Sampling Experiments"
        echo "--"


#         echo "---------------------"
#         echo "Using: $nsamples"

#         printf "\n\n\n"
#         command="$root/$binary
#             --problem=$problem
#             --experiment=tree
#             --ncpus=$ncpus
#             --nsamples=$nsamples
#             --treehight=10
#             --nroots=32"
#         echo $command
#         printf "\n"
#         date
#         printf "\n"
#         time $command

        printf "\n\n\n"
        command="$root/$binary
            --problem=$problem
            --experiment=tree
            --ncpus=$ncpus
            --nsamples=$nsamples
            --treehight=50
            --nroots=32"
        echo $command
        printf "\n"
        date
        printf "\n"
%        time $command



#         printf "\n\n\n"
#         command="$root/$binary
#             --problem=$problem
#             --experiment=tree
#             --ncpus=$ncpus
#             --nsamples=$nsamples
#             --treehight=10000
#             --nroots=1"
#         echo $command
#         printf "\n"
#         date
#         printf "\n"
#         time $command


        printf "\n\n\n"
        command="$root/$binary
            --problem=$problem
            --experiment=async
            --ncpus=$ncpus
            --nsamples=$nsamples"
        echo $command
        printf "\n"
        date
        printf "\n"
%        time $command



        
    done
done



printf "\n\n\n"
command="$root/$binary
            --problem=$problem
            --experiment=async
            --ncpus=16
            --nsamples=100000"
echo $command
printf "\n"
date
printf "\n"
%time $command

