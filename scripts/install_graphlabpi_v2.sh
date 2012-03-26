#!/bin/bash

hg clone https://code.google.com/p/graphlabapi/ 
cd graphlabapi
hg up v2
./configure --bootstrap \
    --force_cmake_install \
    --force_boost_install \
    --eigen \
    --tcmalloc --force_tcmalloc_install \
    --yes | tee mac_configure_log.txt
cd debug
make -j2 | tee ../mac_debug_log.txt

