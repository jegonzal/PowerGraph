#!/bin/bash

create_powerlaw_graph=$PWD/create_powerlaw_graph
convert_graph=$PWD/convert_graph

nverts=10000000

alphas=$(seq 1.8 .2 3)


mkdir datasets
cd datasets

function make_graphs {
    alpha=$1
    dirname=alpha_${alpha}
    echo $dirname
    mkdir $dirname
    cd $dirname
    echo "Making graph"
    graphname=nverts_${nverts}_fanout_alpha_${alpha}
    echo $graphname
    $create_powerlaw_graph --nverts=$nverts --alpha=$alpha \
        --graph ${graphname}.tsv
    echo "Constructing adj"
    $convert_graph ${graphname}.tsv edge2adj ${graphname}.adj
    
    echo "Reversing graph"
    rev_graphname=nverts_${nverts}_fanin_alpha_${alpha}
    echo $graphname   
    $convert_graph ${graphname}.tsv reverse ${rev_graphname}.tsv
    echo "constructing adj"
    $convert_graph ${rev_graphname}.tsv edge2adj ${rev_graphname}.adj  
    
    mkdir ${graphname}_tsv ${graphname}_adj \
        ${rev_graphname}_tsv ${rev_graphname}_adj 

    gzip ${graphname}.tsv 
    gzip ${graphname}.adj 
    gzip ${rev_graphname}.tsv 
    gzip ${rev_graphname}.adj 

    mv ${graphname}.tsv.gz      ${graphname}_tsv/.
    mv ${rev_graphname}.tsv.gz  ${rev_graphname}_tsv/.
    mv ${graphname}.adj.gz      ${graphname}_adj/.
    mv ${rev_graphname}.adj.gz  ${rev_graphname}_adj/.
    cd ..
}


for alpha in $alphas; do
    make_graphs $alpha &
done
