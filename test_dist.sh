#!/bin/sh

#script for auto test graphlab distribution, written by danny bickson
if [ $# -ne 1 ]; then
   echo "Usage: $0 <release number>"
fi
cd /tmp/
rm -fR graphlabapi*
wget http://graphlabapi.googlecode.com/files/graphlabapi_v1_$1.tar.gz
tar xvzf graphlabapi_v1_$1.tar.gz
cd /tmp/graphlabapi
echo "Y" | ./configure --bootstrap
cd release
make -j8


