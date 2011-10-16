#!/bin/bash


hosts=( $( cat $1 ) )
numhosts=${#hosts[@]}
files=( $( ls $2 ) )
numfiles=${#files[@]}
dest_path=$3

echo $hosts
echo $numhosts
echo $files
echo $numfiles
echo $dest_path

for i in $(seq 1 $numfiles) 
do
    let "hostid= $i % $numhosts"
    command="scp ${files[$i]} ${hosts[$hostid]}:$dest_path"
    echo $command
    eval $command
done


