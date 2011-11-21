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
echo $PWD | grep debug > /dev/null
dbgpath=$?
echo $PWD | grep release > /dev/null
relpath=$?
echo $PWD | grep profile > /dev/null
propath=$?

if [ $dbgpath -eq 1 ]; then
  if [ $relpath -eq 1 ]; then
    if [ $propath -eq 1 ]; then
	echo "This test must be run from either ./release/tests/, ./debug/tests/, or ./profile/tests/ in Graphlab root folder"
        echo "Please compile GraphLab first, using the instructions on http://graphlab.org/download.html and try again from the approprite folder"
        exit 1
    fi 
  fi
fi



rm -f $stdoutfname $stderrfname

if [ $# -eq 0 ]; then

echo "Running Standard unit tests"
echo "==========================="
ctest -O testlog.txt
./anytests
./anytests_loader
# delete extra generated files
rm -f dg*

fi

echo 
echo "Running application tests"
echo "========================="
echo "GraphLab collaborative filtering library"
somefailed=0
if [ -f ../demoapps/pmf/pmf ] && [ -f ../demoapps/pmf/itdiff ]; then
  pushd . > /dev/null
  cd ../demoapps/pmf
  echo "---------PMF-------------" >> $stdoutfname
  OUTFILE=smalltest.out
  ./pmf --show_version=true
  if [ $? -eq 2 ]; then
    echo "detected Eigen based pmf"
    OUTFILE=smalltest_eigen.out
  else
    echo "detected it++ based pmf"
  fi
  rm -f smalltest-20-21.out
  echo "********************TEST1************************" >> $stdoutfname
  ./pmf smalltest 0 --scheduler="round_robin(max_iterations=20,block_size=1)" --ncpus=1 --float=true --debug=true >> $stdoutfname 2>& 1 
  
  if ./itdiff smalltest-20-21.out $OUTFILE ; then
    echo "PASS TEST 1 (Alternating least sqaures)"
  else
     somefailed=1
    echo "FAIL: Output differs!"
  fi
  echo "********************TEST2************************" >> $stdoutfname
  ./pmf --unittest 1 --ncpus=1 --debug=true >> $stdoutfname 2>& 1
  if [ $? -eq 0 ]; then
     echo "PASS TEST 2 (Alternating least squares)"
  else
     somefailed=1
     echo "FAIL --unittest=1"
  fi
  echo "********************TEST3************************" >> $stdoutfname
  ./pmf --unittest 71 --ncpus=1 --debug=true >> $stdoutfname 2>& 1
  if [ $? -eq 0 ]; then
     echo "PASS TEST 3 (Lanczos)"
  else
     somefailed=1
     echo "FAIL --unittest=71 (Lanczos)"
  fi
  echo "********************TEST4************************" >> $stdoutfname
  ./pmf --unittest 91 --ncpus=1 --debug=true >> $stdoutfname 2>& 1
  if [ $? -eq 0 ]; then
     echo "PASS TEST 4 (Weighted ALS)"
  else
     somefailed=1
     echo "FAIL --unittest=91 (weighted alternating least squares)"
  fi
  echo "********************TEST5************************" >> $stdoutfname
 ./pmf --unittest 101 --ncpus=1 >> $stdoutfname 2>& 1 
  if [ $? -eq 0 ]; then
     echo "PASS TEST 5 (CoSaMP)"
  else
     echo "FAIL --unittest=101 (CoSaMP)"
     somefailed=1
  fi
  echo "********************TEST6************************" >> $stdoutfname
 ./pmf --unittest 131  >> $stdoutfname 2>& 1 
  if [ $? -eq 0 ]; then
     echo "PASS TEST 6 (SVD)"
  else
     echo "FAIL --unittest=131 (SVD)"
     somefailed=1
  fi
  popd > /dev/null

else
  echo "PMF not found. "
fi
echo



echo "GraphLab clustring library"
if [ -f ../demoapps/clustering/glcluster ]; then
  pushd . > /dev/null
  cd ../demoapps/clustering
  echo "---------CLUSTERING-------------" >> $stdoutfname
  echo "---------CLUSTERING-------------" >> $stderrfname
  echo "********************TEST1************************" >> $stdoutfname
  ./glcluster --unittest 1  $stdoutfname 2>& 1
  if [ $? -eq 0 ]; then
     echo "PASS TEST 1 (Math functions)"
  else
     somefailed=1
     echo "FAIL --unittest=1 (Math functions)"
  fi
  echo "********************TEST2************************" >> $stdoutfname
  ./glcluster --unittest 2 >> $stdoutfname 2>& 1
  if [ $? -eq 0 ]; then
     echo "PASS TEST 2 (Distance functions)"
  else
     somefailed=1
     echo "FAIL --unittest=2 (Distance functions)"
  fi
  echo "********************TEST3************************" >> $stdoutfname
  ./glcluster --unittest 4 >> $stdoutfname 2>& 1
  if [ $? -eq 0 ]; then
     echo "PASS TEST 3 (Floating point math functions)"
  else
     somefailed=1
     echo "FAIL --unittest=3 (Floating point math functions)"
  fi
else
  echo "Clustering library not found. "
fi
 
  popd  > /dev/null

  if [ $somefailed == 1 ]; then
     echo "**** FAILURE LOG **************" >> $stdoutfname
     cat $stderrfname >> $stdoutfname
     echo "**** CONFIGURE.DEPS **************" >> $stdoutfname
     cat ../../configure.deps >> $stdoutfname
     echo "**** CONFIG.LOG **************" >> $stdoutfname
     cat ../../config.log >> $stdoutfname
     echo "**** SYSTEM STATS **************" >> $stdoutfname
     echo `date` >> $stdoutfname
     echo `uname -a` >> $stdoutfname
     echo `echo $USER` >> $stdoutfname
     echo "Some of the tests failed".
     echo "Please email stdout.log to danny.bickson@gmail.com"
     echo "Thanks for helping improve GraphLab!"
  fi




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
