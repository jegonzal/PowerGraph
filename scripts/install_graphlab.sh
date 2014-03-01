#!/bin/bash

git clone https://code.google.com/p/graphlabapi/ 
cd graphlab
./configure | tee install_configure_log.txt
cd release
make -j2 | tee ../v2_debug_log.txt

