#!/bin/bash
IN_BOOTSTRAP=1
mkdir -p deps
cd deps/
source ../cmakecheck.sh
source ../boostcheck.sh
cd ..
