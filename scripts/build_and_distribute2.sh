#!/bin/bash
cd ~/graphlabapi/release/apps/dist2_apps/webgraph_bp
make -j8 webgraph_bp

cp ~/graphlabapi/release/apps/dist2_apps/webgraph_bp/webgraph_bp /mnt/bin/

cd /mnt
mpiexec -n 32 -f ~/hosts32 -wdir /mnt /mnt/bin/dirdist -s bin -f 3
