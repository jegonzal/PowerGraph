#!/bin/bash


for trial in `seq 1 10`; do
    printf "\n\n\n"
    command="../pgibbs_jt_alchemy
            ./elidan1.out
            --treesize=1000 
            --ncpus=8
            --runtime={1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,150,200,250,300,350,400,500}
            --treewidth=4"
    echo $command
    printf "\n"
    date
    printf "\n"
    time eval $command
    
    command="../pgibbs_jt_alchemy
            ./elidan1.out
            --treesize=1000 
            --ncpus=8
            --runtime={1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,150,200,250,300,350,400,500}
            --treewidth=4
            --priorities=true"
    echo $command
    printf "\n"
    date
    printf "\n"
    time eval $command
done