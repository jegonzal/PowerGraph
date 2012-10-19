#!/bin/bash
if [ ! -d release ]; then
  echo "Run from the graphlab root folder after ./configure"
else
  cd release
  make external_dependencies
  cd ..
  rm -f deps/local/lib/libboost*.so deps/local/lib/libhdfs*.so deps/local/lib/libtcmalloc*.so deps/local/lib/libevent*.so
  rm -f deps/local/lib/libboost*.dylib deps/local/lib/libhdfs*.dylib deps/local/lib/libtcmalloc*.dylib deps/local/lib/libevent*.dylib
  cd release
  make -j2
fi
