#!/bin/bash


hosts=( $( cat "hosts" ) )
numhosts=${#hosts[@]}
fileprefix=altavista
dest_path=/mnt/webgraph

echo $hosts
echo $numhosts
echo $files
echo $numfiles
echo $dest_path

for i in $(seq 0 511)
do
    let "hostid= $i % $numhosts"
    filename=$fileprefix.$i.dump
    command="scp /mnt/external/webgraphs/full_graph/$filename ${hosts[$hostid]}:$dest_path/. &"
    echo $command
    eval $command
done


