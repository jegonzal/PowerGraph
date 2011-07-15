#!/bin/bash

pushd .
cd $DEPS_DIR

function test_for_avro {
  # if force avro install is set, don't probe. assume failure
  if [ ! -z $force_avro_install ] ; then
    return
  fi
  # list all the avro files we use
  grep -h -r "include <avro/*" ../src | uniq > avro_tester.cpp
  # make a main which does nothing
  echo -e "\nint main(int argc, char** argv){ return 0; }\n" >> avro_tester.cpp
  # try to compile it. Use AVRO_ROOT if available
  rm -f a.out
  if [ -z $AVRO_ROOT ] ; then
    echo "AVRO_ROOT not defined. Probing in usual directories..."
    g++ -L/usr/local/lib -I/usr/local/include -lavrocpp avro_tester.cpp > /dev/null 2> /dev/null
  else
    echo "Probing in $AVRO_ROOT"
    g++ -L$AVRO_ROOT/lib -I$AVRO_ROOT/include -lavrocpp avro_tester.cpp > /dev/null 2> /dev/null
  fi
  if [ -f a.out ] ; then
    avrofound=1
  fi
}


function download_file_with_forward {
  # detect wget
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


if [ -z $IN_BOOTSTRAP ]; then
  echo "This script should not be run directly."
else
  echo "Probing for avro..."
  installprefix=$PWD
  test_for_avro
  # test if there is an existing installation
  if [ ! -z $force_avro_install ] || [ -z $avrofound ]; then
    AVRO_ROOT=$installprefix
    test_for_avro
    if [ ! -z $avrofound ]; then
      echo "Existing installation of avro found in $AVRO_ROOT"
    fi
  fi
  if [ ! -z $force_avro_install ] || [ -z $avrofound ]; then

    echo " ==================== Avro Not Found ! ==================== "
    echo
    echo "Avro libraries were not detected. Either your system does not have"
    echo "Avro, or AVRO_ROOT was not specified correctly."
    echo "The script will now proceed to download boost and install it locally"
    echo "in the graphlabapi/deps directory."
    echo
    echo "Press Enter to proceed,"
    echo "or Ctrl-C to break and install AVRO yourself, or try again with AVRO_ROOT"
    if [ -z $always_yes ] ; then
      read
    fi
    echo "Download Boost 1.5.1 from UOregon Mirror..."
    #download boost
    rm avro-cpp-1.5.1.tar.gz
    download_file_with_forward \
        "http://mirror.uoregon.edu/apache//avro/avro-1.5.1/cpp/avro-cpp-1.5.1.tar.gz" \
        "avro-cpp-1.5.1.tar.gz"
    # If there is already an avro 1.5.1 here. delete it
    echo "Cleaning up prior installations of boost..."
    if [ -d "avro_1_5_1" ] ; then
      rm -rf avro_1_5_1
    fi
    if [ -d "include/avro" ] ; then
      rm -rf include/avro
    fi
    if [ -d "lib" ] ; then
      rm -f lib/libavrocpp*
    fi
    set -e
    tar -xjvf avro-cpp-1.5.1.tar.gz
    # cd into avro directory
    cd avro-cpp-1.5.1
    # configure avro cmake
    echo "Configuring Avro..."
    $CMAKE "-DCMAKE_INSTALL_PREFIX:STRING=$installprefix"
    echo "Building Avro"
    make | tee avrolog.txt
    echo "Installing Avro"
    make install
    cd ..
    set +e
    # set the AVRO_ROOT
    AVRO_ROOT=$installprefix
  else
    echo "Avro Found!"
  fi
fi

popd



