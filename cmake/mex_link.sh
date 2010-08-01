#!/bin/sh

OUTPUT=$1
STATIC_LIB_NAME=$OUTPUT
MOVE_LOCATION=$2
BASEDIR=$3
LINKFILES=$4
LINKER_FLAGS= -shared -Wl,--version-script,/afs/cs.cmu.edu/misc/matlab/amd64_f7/7.9/lib/matlab7/extern/lib/glnxa64/mexFunction.map -Wl,--no-undefined

echo $OUTPUT
echo $STATIC_LIB_NAME

 

echo 'Running godawful mex linking hack...'
# mex_stub.o compiled with:
g++ -g -Wall -fPIC -ansi -D_GNU_SOURCE -fPIC -fno-omit-frame-pointer  -pthread  -DMATLAB_MEX_FILE -lmx -lmex -lmat -lm -I/afs/cs.cmu.edu/local/matlab/amd64_f7/7.9/lib/matlab7/extern/include -c ${BASEDIR}/cmake/Mex_stub.cpp -o mex_stub.o
mex -g -cxx CC='gcc' CXX='g++' LD='g++' -L./ -lglib-2.0 -l$STATIC_LIB_NAME $LINKER_FLAGS -output $OUTPUT mex_stub.o $LINKFILES
#mv $OUTPUT.mexa64 $MOVE_LOCATION

#-lpthread -lgthread-2.0 -lrt -DMX_COMPAT_32