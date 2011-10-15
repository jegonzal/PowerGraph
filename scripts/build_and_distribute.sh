#!/bin/bash

# hostfile="~/hosts"
hostfile="~/hosts32"

cd ~/graphlabapi/release/apps/dist2_apps/webgraph_bp
make -j8
cp * ~/bin

for i in `cat ~/hosts32`; do
    echo "connecting to $i";
    rsync -avz ~/bin $i:~/.;
done
