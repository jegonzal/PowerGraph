#!/bin/bash
# create an empty configure.deps
echo "CMAKE='cmake'" > configure.deps
IN_BOOTSTRAP=1
mkdir -p deps
cd deps/
source ../cmakecheck.sh
source ../boostcheck.sh
cd ..