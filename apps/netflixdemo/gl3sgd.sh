#!/bin/sh
cmd="mpiexec -np 12 -hostfile $HOME/machines ../../release/apps/netflixdemo/gl3sgd_netflix_demo --matrix=$HOME/data/netflix/split --movielist=$HOME/data/netflix/moviename.txt --iterations=3 --D=20"
echo $cmd
$cmd
