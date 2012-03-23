#!/bin/bash

create_powerlaw_graph=$PWD/create_powerlaw_graph
convert_graph=$PWD/convert_graph

nverts=10000000
fanout=10
alphas=$(seq 1 .2 3)


mkdir datasets
cd datasets


for alpha in $alphas; do
    dirname=alpha_${alpha}
    echo $dirname
    mkdir $dirname
    cd $dirname
    echo "Making graph"
    graphname=nverts_${nverts}_fanout_${fanout}_alpha_${alpha}
    echo $graphname
    $create_powerlaw_graph --nverts=$nverts --fanout=$fanout --alpha=$alpha \
        --graph ${graphname}.tsv
    echo "Constructing adj"
    $convert_graph ${graphname}.tsv edge2adj ${graphname}.adj
    
    echo "Reversing graph"
    rev_graphname=nverts_${nverts}_fanin_${fanout}_alpha_${alpha}
    echo $graphname   
    $convert_graph ${graphname}.tsv reverse ${rev_graphname}.tsv
    echo "constructing adj"
    $convert_graph ${rev_graphname}.tsv edge2adj ${rev_graphname}.adj  
    cd ..
done
