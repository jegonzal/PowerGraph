#!/bin/bash
echo "Running Standard unit tests"
echo "==========================="
#ctest -O testlog.txt
anytests
anytests_loader

echo 
echo "Running application tests"
echo "========================="
echo "PMF..."
cd ../demoapps/pmf
./pmf  smalltest 0 --scheduler="round_robin(max_iterations=10)" --float=true --ncpus=1 > /dev/null 2> /dev/null
if ./itdiff smalltest20.out smalltest.out ; then
  echo "PASS"
else
  echo "FAIL: Output differs!"
  exit 1
fi

echo

cd ../demo
echo "Demo..."
./demo > /dev/null 2> /dev/null
if [ $? -eq 0 ]; then
  echo "PASS"
else
  echo "FAIL. Demo App returned with failure"
  exit 1
fi