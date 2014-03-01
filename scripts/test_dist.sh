#!/bin/sh

#script for auto test graphlab distribution, written by danny bickson
if [ $# -ne 1 ]; then
   echo "Usage: $0 <release number>"
fi
rm -fR /tmp/graphlabapi*
cp dist/graphlabapi_v1_$1.tar.gz /tmp/
cd /tmp/

tar xvzf graphlabapi_v1_$1.tar.gz
cd /tmp/graphlabapi
./configure --bootstrap --yes
cd release
make -j8

cd tests
./runtests.sh
