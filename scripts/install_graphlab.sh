#!/bin/bash

hg clone https://code.google.com/p/graphlabapi/ 
cd graphlabapi
hg up v2
./configure | tee v2_install_configure_log.txt
cd debug
make -j2 | tee ../v2_debug_log.txt

