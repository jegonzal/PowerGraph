#!/bin/bash

function quit_if_bad_retvalue {
  if [ $? -eq 0 ]; then
    echo "PASS"
  else
    echo "FAIL. Program returned with failure"
    exit 1
  fi
}

function test_rpc_prog {
  echo "Testing $1 ..."
  mpiexec -n 2 -host $localhostname ./$1 > /dev/null 2> /dev/null
  if [ $? -ne 0 ]; then
    echo "FAIL. Program returned with failure"
    exit 1
  fi
  str="mpiexec -n 2 -host $localhostname ./$1 2> /dev/null | grep \"$2\""
  echo $str
  e=`eval $str`
  if [ -z "$e" ] ; then
    echo "Expected program output not obtained"
    exit 1
  fi
}


echo "Running Standard unit tests"
echo "==========================="
ctest -O testlog.txt
./anytests
./anytests_loader
# delete extra generated files
rm -f dg*

pushd . > /dev/null

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
quit_if_bad_retvalue

popd > /dev/null

echo "RPC Tests"
echo "========="
echo "Testing for availability of an MPI daemon"
localhostname=`hostname`
mpdtrace
if [ $? -eq 0 ]; then
  echo "MPI available"
else
  echo "MPI not available. Distributed/RPC tests not running."
  exit 1
fi



test_rpc_prog rpc_example1 "5 plus 1 is : 6\\|11 plus 1 is : 12"
test_rpc_prog rpc_example2 "hello world!\\|1, 2, 1,"
test_rpc_prog rpc_example3 "1.a = 10\\|10.b = 0\\|string = hello world!"
test_rpc_prog rpc_example4 "1.a = 10\\|10.b = 0\\|string = hello world!"
test_rpc_prog rpc_example5 "1 + 2.000000 = three"
test_rpc_prog rpc_example6 "10\\|15\\|hello world\\|10.5\\|10"
test_rpc_prog rpc_example7 "set from 1\\|set from 1\\|set from 0\\|set from 0\\|set from 1\\|set from 1\\|set from 0\\|set from 0"

echo "Testing Distributed disk graph construction..."
mpiexec -n 2 -host $localhostname ./distributed_dg_construction_test > /dev/null 2> /dev/null
quit_if_bad_retvalue
rm dg*