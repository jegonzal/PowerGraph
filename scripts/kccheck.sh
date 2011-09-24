#!/bin/bash

pushd .
cd $DEPS_DIR

function test_for_kc {
  # if force boost install is set, don't probe. assume failure
  if [ ! -z $force_kc_install ] ; then
    return
  fi
  # make a main which does nothing
  echo -e "#include <kchashdb.h>\n" > kc_tester.cpp
  echo -e "\nint main(int argc, char** argv){ return 0; }\n" >> kc_tester.cpp
  # try to compile it. Use KC_ROOT if available
  rm -f a.out
  if [ -z $KC_ROOT ] ; then
    echo "KC_ROOT not defined. Probing in usual directories..."
    g++ -L/usr/local/lib -I/usr/local/include -lkyotocabinet kc_tester.cpp > /dev/null 2> /dev/null
  else
    echo "Probing in $KC_ROOT"
    g++ -L$KC_ROOT/lib -I$KC_ROOT/include -lkyotocabinet kc_tester.cpp > /dev/null 2> /dev/null
  fi
  if [ -f a.out ] ; then
    kcfound=1
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
  echo "Probing for Kyoto Cabinet..."
  installprefix=$PWD
  test_for_kc
  # test if there is an existing installation
  if [ ! -z $force_kc_install ] || [ -z $kcfound ]; then
    KC_ROOT=$installprefix
    test_for_kc
    if [ ! -z $kcfound ]; then
      echo "Existing installation of kc found in $KC_ROOT"
    fi
  fi
  if [ ! -z $force_kc_install ] || [ -z $kcfound ]; then

    echo " ==================== Kyoto Cabinet Not Found ! ==================== "
    echo
    echo "Kyoto Cabinet was not detected. Either your system does not have"
    echo "Kyoto Cabinet, or KC_ROOT was not specified correctly."
    echo "The script will now proceed to download Kyoto Cabinet and install it"
    echo "locally in the graphlabapi/deps directory."
    echo
    echo "Press Enter to proceed,"
    echo "or Ctrl-C to break and install Kyoto Cabinet yourself, or try again "
    echo "with KC_ROOT"
    if [ -z $always_yes ] ; then
      read
    fi
    echo "Download Kyoto Cabinet 1.2.70 from http://fallabs.com/  ..."
    #download kyoto cabinet
    rm kyotocabinet.tar.gz
    download_file_with_forward http://fallabs.com/kyotocabinet/pkg/kyotocabinet-1.2.70.tar.gz kyotocabinet.tar.gz
    # If there is already a kyoto cabinet here. delete it
    echo "Cleaning up prior installations of Kyoto Cabinet..."
    if [ -d "kyotocabinet-1.2.70" ] ; then
      rm -rf kyotocabinet-1.2.70
    fi
    if [ -d "include" ] ; then
      rm -rf include/kc*.h
    fi
    if [ -d "lib" ] ; then
      rm -f lib/libkyotocabinet*
    fi
    set -e
    tar -xzvf kyotocabinet.tar.gz
    # cd into kyoto cabinet directory
    cd kyotocabinet-1.2.70
    # build kyoto cabinet
    echo "Configuring Kyoto Cabinet..."
    ./configure --prefix=$installprefix
    echo "Compiling Kyoto Cabinet..."
    make -j2
    
    echo "Installing Kyoto Cabinet..."
    make install
    cd ..
    set +e
    # set KC_ROOT
    KC_ROOT=$installprefix
  else
    echo "Kyoto Cabinet Found!"
  fi
fi

popd

