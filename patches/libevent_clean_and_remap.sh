#!/bin/bash
if [ $# -ne 1 ]
then
  echo "Missing argument: directory where libevent libraries reside"
  exit 1
fi

CURDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $1
rm -f libevent*.so
objcopy --redefine-syms=$CURDIR/libevent_remap_file.txt libevent_pthreads.a
objcopy --redefine-syms=$CURDIR/libevent_remap_file.txt libevent.a

