#!/bin/bash
cd /mnt/graphlabapi/release/apps/dist2_apps/webgraph_bp
make -j8 webgraph_bp
cp /mnt/graphlabapi/release/apps/dist2_apps/webgraph_bp/webgraph_bp /mnt/bin
cd /mnt
mpiexec -n 16 -f ~/hosts -wdir /mnt /mnt/bin/dirdist -s bin -f 3
mpiexec -n 16 -f ~/hosts chmod +x /mnt/bin/*



