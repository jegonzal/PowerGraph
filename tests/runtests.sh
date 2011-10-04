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
  echo "---------$1-------------" >> $stdoutfname
  echo "---------$1-------------" >> $stderrfname 
  mpiexec -n 2 -host $localhostname ./$1  >> $stdoutfname 2>> $stderrfname
  if [ $? -ne 0 ]; then
    echo "FAIL. Program returned with failure"
    exit 1
  fi
  str="mpiexec -n 2 -host $localhostname ./$1 2> /dev/null | grep \"$2\""
  #echo $str
  e=`eval $str`
  if [ -z "$e" ] ; then
    echo "Expected program output not obtained"
    exit 1
  fi
}

stdoutfname=$PWD/stdout.log
stderrfname=$PWD/stderr.log
rm -f $stdoutfname $stderrfname

echo "Running Standard unit tests"
echo "==========================="
ctest -O testlog.txt
./anytests
./anytests_loader
# delete extra generated files
rm -f dg*


echo 
echo "Running application tests"
echo "========================="
echo "PMF..."
if [ -f ../demoapps/pmf/pmf ] && [ -f ../demoapps/pmf/itdiff ]; then
  pushd . > /dev/null
  cd ../demoapps/pmf
  echo "---------PMF-------------" >> $stdoutfname
  echo "---------PMF-------------" >> $stderrfname
  ./pmf  smalltest 0 --scheduler="round_robin(max_iterations=10)" --float=true --ncpus=1  >> $stdoutfname 2>> $stderrfname
  if ./itdiff smalltest-20-21.out smalltest.out ; then
    echo "PASS TEST 1 (Alternating least sqaures)"
  else
    echo "FAIL: Output differs!"
    exit 1
  fi

  ./pmf --unittest 1 --ncpus=1 >> $stdoutfname 2>> $stderrfname 
  if [ $? -eq 0 ]; then
     echo "PASS TEST 2 (alternating least squares)"
  else
     echo "FAIL --unittest=1"
     exit 1
  fi
 ./pmf --unittest 91 --ncpus=1 --debug=true >> $stdoutfname 2>> $stderrfname 
  if [ $? -eq 0 ]; then
     echo "PASS TEST 3"
  else
     echo "FAIL --unittest=91 (weighted alternating least squares)"
     exit 1
  fi
 ./pmf --unittest 101 --ncpus=1 >> $stdoutfname 2>> $stderrfname 
  if [ $? -eq 0 ]; then
     echo "PASS TEST 4"
  else
     echo "FAIL --unittest=101 (CoSaMP)"
     exit 1
  fi
   
  popd  > /dev/null
else
  echo "PMF not found. "
fi
echo

if [ -f ../demoapps/demo/demo ]; then
  pushd . > /dev/null
  cd ../demoapps/demo
  echo "Demo..."
  echo "---------demo-------------" >> $stdoutfname
  echo "---------demo-------------" >> $stderrfname
  
  ./demo  >> $stdoutfname 2>> $stderrfname
  quit_if_bad_retvalue
  popd > /dev/null
else
  echo "demo not found. "
fi

echo
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

echo
echo "Distributed GraphLab Tests"
echo "=========================="

echo "Testing Distributed disk graph construction..."
echo "---------distributed_dg_construction_test-------------" >> $stdoutfname
echo "---------distributed_dg_construction_test-------------" >> $stderrfname 
mpiexec -n 2 -host $localhostname ./distributed_dg_construction_test >> $stdoutfname 2>> $stderrfname
quit_if_bad_retvalue
rm -f dg*

echo "Testing Distributed Graph ..."
echo "---------distributed_graph_test-------------" >> $stdoutfname
echo "---------distributed_graph_test-------------" >> $stderrfname 
./distributed_graph_test -g
mpiexec -n 2 -host $localhostname ./distributed_graph_test -b >> $stdoutfname 2>> $stderrfname
quit_if_bad_retvalue
rm -f dg*
