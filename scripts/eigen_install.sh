#!/bin/bash


_Box () {
    str="$@"
    len=$((${#str}+4))
    for (( i=1; i<=$len; i++ )); do echo -n '-'; done;
    echo; echo "| "$str" |";
    for (( i=1; i<=$len; i++ )); do echo -n '-'; done;
    echo
}

function test_for_eigen {
   echo "force install"
   #if [ ! -d $1/deps/eigen-eigen-3.0.2/ ]; then
   #   eigenfound=1;
   #   echo "Eigen was found at: $1/deps/eigen-eigen-3.0.2/"
   #fi
}


function download_file_with_forward {
  echo "Downloading $2 from $1 ..."
  wgetcmd=`which wget`
  if [ -z $wgetcmd ] ; then
    curlcmd=`which curl`
    if [ -z $curlcmd ] ; then
      echo "Unable to find either curl or wget! Cannot proceed with automatic install."
      exit 1
    fi
    curl -L $1 -o $2
  else
    wget $1 -O $2
  fi
}


function install_from_source {
  if  [ $eigeninstalled -eq 1 ]; then
    return
  fi
  
  echo " ============= Source installation of Eigen =============== "
  # detect wget
  echo "The script will now proceed to download Eigen and try to set it up"
  echo
  echo "Press Enter to proceed,"
  echo "or Ctrl-C to break and install Eigen"
  read
  echo "Downloading Eigen 3.0.2 from www.tuxfamily.com... into `pwd`"
  download_file_with_forward http://bitbucket.org/eigen/eigen/get/3.0.2.tar.gz $1/deps/3.0.2.tar.gz
  # unpack
  set -e
  cd $1/deps/
  tar -xzvf 3.0.2.tar.gz
  set +e

  echo "Eigen setup success!"
}

pushd .
curdir=$1
mkdir -p $curdir/$DEPS_DIR
cd $curdir/$DEPS_DIR
echo "Detecting Eigen... (enered from dir: $curdir).."
pwd
echo
test_for_eigen $1

if [ -z $eigenfound ] ; then
  echo " ==================== Eigen Not Found ! ===================== "
  OS=`uname`
  eigeninstalled=0
  install_from_source $1
fi
popd
