#!/bin/sh
mpiexec -np 16 -hostfile $HOME/machines killall als_netflix_demo 
mpiexec -np 16 -hostfile $HOME/machines killall gl3als_netflix_demo 
