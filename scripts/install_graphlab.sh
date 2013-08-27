#!/bin/bash

git clone https://github.com/graphlab-code/graphlab.git
cd graphlab
./configure | tee install_configure_log.txt
cd release
make -j2 | tee ../v2_build_log.txt

