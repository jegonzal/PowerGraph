#!/bin/sh
mpiexec -np 16 -hostfile $HOME/machines killall als_netflix_demo 
