#!/bin/sh
cmd="mpiexec -np 12 -hostfile $HOME/machines env CLASSPATH=`hadoop classpath` ../../release/apps/netflixdemo/als_netflix_demo --matrix=$HOME/data/netflix/split --movielist=$HOME/data/netflix/moviename.txt --engine=sync --max_iter=1 --D=20 --maxval=5 --minval=0 --lambda=1 --regnormal=0"
#--implicitratingweight=1 --implicitratingvalue=3 --implicitratingpercentage=0.01 --users=480189 --items=17770 --implicitratingtype=1"
echo $cmd
$cmd
