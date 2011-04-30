#!/bin/bash
echo "Running Standard unit tests"
echo "==========================="
ctest -O testlog.txt
anytests
anytests_loader

echo 
echo "Running application tests"
echo "========================="
cd ../demoapps/pmf
./pmf  smalltest 0 --scheduler="round_robin(max_iterations=10)" --float=true --ncpus=1 
if diff smalltest20.out smalltest.out >/dev/null ; then
  echo "PASS"
else
  echo "FAIL: Output differs!"
fi

