#!/bin/bash


hosts=( $( cat $1 ) )
numhosts=${#hosts[@]}
fileprefix=atom
dest_path=/mnt/webgraph

echo $hosts
echo $numhosts
echo $files
echo $numfiles
echo $dest_path

for i in $(seq 0 511)
do
    let "hostid= $i % $numhosts"
    filename=$fileprefix$i
    command="scp $filename ${hosts[$hostid]}:$dest_path/."
    echo $command
#    eval $command
done


