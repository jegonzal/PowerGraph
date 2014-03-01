#!/bin/bash

if [ ! -z "$PMI_RANK" ]; then
RANK=$PMI_RANK
elif [ ! -z "$OMPI_COMM_WORLD_RANK" ]; then
RANK=$OMPI_COMM_WORLD_RANK
else
echo "Unable to figure out MPI Rank!"
exit 1
fi
#echo $RANK
$* 2>&1 | tee out.$RANK 
