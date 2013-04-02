#!/bin/sh
cmd="mpiexec -np 12 -hostfile $HOME/machines ../../release/apps/netflixdemo/als_netflix_demo --matrix=$HOME/data/netflix/split --movielist=$HOME/data/netflix/moviename.txt --engine=sync --max_iter=3 --D=20"
echo $cmd
$cmd
