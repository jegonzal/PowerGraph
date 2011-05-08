#!/bin/bash
function test_for_boost {
  # if force boost install is set, don't probe. assume failure
  if [ ! -z $force_boost_install ] ; then
    return
  fi
  # list all the boost files we use
  grep -h -r "include <boost/*" ../src | uniq > boost_tester.cpp
  # make a main which does nothing
  echo -e "\nint main(int argc, char** argv){ return 0; }\n" >> boost_tester.cpp
  # try to compile it. Use BOOST_ROOT if available
  rm -f a.out
  if [ -z $BOOST_ROOT ] ; then
    echo "BOOST_ROOT not defined. Probing in usual directories..."
    g++ -L/usr/local/lib -I/usr/local/include -lboost_program_options -lboost_filesystem -lboost_system boost_tester.cpp
  else
    echo "Probing in $BOOST_ROOT"
    g++ -L$BOOST_ROOT/lib -I$BOOST_ROOT/include -lboost_program_options -lboost_filesystem -lboost_system boost_tester.cpp
  fi
  if [ -f a.out ] ; then
    boostfound=1
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
  echo "Probing for boost..."
  installprefix=$PWD
  test_for_boost
  # test if there is an existing installation
  if [ ! -z $force_boost_install ] || [ -z $boostfound ]; then
    BOOST_ROOT=$installprefix
    test_for_boost
    if [ ! -z $boostfound ]; then
      echo "Existing installation of boost found in $BOOST_ROOT"
    fi
  fi
  if [ ! -z $force_boost_install ] || [ -z $boostfound ]; then

    echo " ==================== Boost Not Found ! ==================== "
    echo
    echo "Boost libraries were not detected. Either your system does not have"
    echo "Boost, or BOOST_ROOT was not specified correctly."
    echo "The script will now proceed to download boost and install it locally"
    echo "in the graphlabapi/deps directory."
    echo
    echo "Press Enter to proceed,"
    echo "or Ctrl-C to break and install Boost yourself, or try again with BOOST_ROOT"
    read
    echo "Download Boost 1.46.1 from SourceForge..."
    #download boost
    rm boost.tar.bz2
    download_file_with_forward http://sourceforge.net/projects/boost/files/boost/1.46.1/boost_1_46_1.tar.bz2/download boost.tar.bz2
    # If there is already a boost 1_46_1 here. delete it
    echo "Cleaning up prior installations of boost..."
    if [ -d "boost_1_46_1" ] ; then
      rm -rf boost_1_46_1
    fi
    if [ -d "include/boost" ] ; then
      rm -rf include/boost
    fi
    if [ -d "lib" ] ; then
      rm -f lib/libboost*
    fi
    tar -xjvf boost.tar.bz2
    # cd into boost directory
    cd boost_1_46_1
    # build boost
    echo "Bootstrapping Boost..."
    ./bootstrap.sh --with-libraries="filesystem,program_options,system" --prefix=$installprefix
    echo "Compiling Boost..."
    ./bjam --threading=multi --link=static --variant=release | tee bjamlog.txt
    if ! grep -q "were successfully built" bjamlog.txt ; then
      echo "Boost failed to build for unknown reasons. Installation cannot proceed."
      exit 1
    fi
    echo "Installing Boost...This could take a little while..."
    ./bjam install
    cd ..
    # set BOOST_ROOT
    BOOST_ROOT=$installprefix
  else
    echo "Boost Found!"
  fi
fi
