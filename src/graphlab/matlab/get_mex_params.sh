#!/bin/bash
#
# This file should not be called from anywhere else but the 
# generate_mex_makefile.m matlab function

if [ $# -gt 0 ]
then
  Arch=$ARCH
  MAPFILE="mexFunction.map"
  source $MATLAB/bin/mexopts.sh
  extra_includes="-I\"${MATLAB}/extern/include\""
  if [ -d $MATLAB/simulink ]; then
    extra_includes="${extra_includes} -I\"$MATLAB/simulink/include\""
 fi
  echo "MEXLDFLAGS = $LDFLAGS $CXXLIBS" > $1
  echo "MEXCXXFLAGS = -DMATLAB_MEX_FILE $extra_includes $CXXFLAGS" >> $1
fi

