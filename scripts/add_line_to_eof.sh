#!/bin/bash
for f in `find src \( -name "*.cpp" -or -name "*.hpp" \)`; do
lastline=`tail -n 1 $f`
len=$((${#lastline}))
if  [ $len -ne 0 ]; then
    echo $f
    echo -e "" >> $f
fi
done
