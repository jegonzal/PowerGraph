#!/bin/bash


for time in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 75 100 200 300; do
    for trial in `seq 1 10`; do
        printf "\n\n\n"
        command="../pgibbs_jt_alchemy
            ./original/image.alchemy 
            --treesize=1000 
            --ncpus=16
            --runtime=$time
            --treewidth=3"
        echo $command
        printf "\n"
        date
        printf "\n"
        time eval $command

        command="../pgibbs_jt_alchemy
            ./original/image.alchemy 
            --treesize=1000 
            --ncpus=16
            --runtime=$time
            --treewidth=3
            --priorities=true"
        echo $command
        printf "\n"
        date
        printf "\n"
        time eval $command

        
    done
done