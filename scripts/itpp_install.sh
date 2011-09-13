#!/bin/bash

function find_itpp_lib {
  # if force itpp install is set, don't probe. assume failure
  if [ ! -z $force_itpp_install ] ; then
    return
  fi
  
  # look for itppconfig
  itppconfig_pos=`which itpp-config`
  # if not found and itpp root is defined. try the itpp_root/bin
  if [ ! -z $ITPP_ROOT ] && [ -z $itppconfig_pos ] && [ -f $ITPP_ROOT/bin/itpp-config ]; then
    itppconfig_pos=$ITPP_ROOT/bin/itpp-config
  fi
  
  if [ -z $itppconfig_pos ] ; then
#    echo "itpp-config not found."
    return
  fi
  # itpp-config exists. Probe for its version and location
  ITPPVERSION=`$itppconfig_pos --version`
  ITPPLIB=`$itppconfig_pos --libs --static`
  if [ ! -z $ITPPVERSION ] ; then
    echo "itpp $ITPPVERSION detected at at $ITPPLIB"
  else
#    echo "Failed to obtain itpp version"
    return
  fi
}

_Box () {
    str="$@"
    len=$((${#str}+4))
    for (( i=1; i<=$len; i++ )); do echo -n '-'; done;
    echo; echo "| "$str" |";
    for (( i=1; i<=$len; i++ )); do echo -n '-'; done;
    echo
}

function run_as_sudo {

  echo "To complete the installation process, we will require root access" 
  echo "to execute the following command:"
  echo
  _Box $1
  echo 
  echo "sudo will be used to execute the command and prompt for a password."
  echo
  echo "Press Enter to continue, "
  echo "or press Ctrl-C to break and run the above command yourself."
  read

  sudo -k
  sudo $1
}


function test_for_itpp {

  find_itpp_lib
  echo
#  echo "Trying to link against itpp..."
  # try to link against it
  echo -e "#include <itpp/itbase.h>\n" > itpp_tester.cpp
  echo -e "\nint main(int argc, char** argv){ return 0; }\n" >> itpp_tester.cpp

  rm -f a.out
  
  g++ $ITPPLIB itpp_tester.cpp > /dev/null 2> /dev/null
  
  if [ -f a.out ] ; then
    echo "Probe successful. ITPP should be functional"
    itppfound=1
  else
    g++ $ITPPLIB itpp_tester.cpp -llapack -lblas > /dev/null 2> /dev/null
    if [ -f a.out ] ; then
    echo "Probe successful. ITPP should be functional"
      itppfound=1
    else
    #echo "ITPP not found."
    if [ ! -z $ITPPVERSION ]; then
      echo "There is a problem with your itpp installation."
      echo "itpp-config was found, but we are unable to link against itpp."
    fi
    fi
  fi
}


function try_install_from_apt_get {
  if  [ $itppinstalled -eq 1 ]; then
    return
  fi
  
  # ok search for apt-get
  aptget=`which apt-get`
  if [ -z $aptget ] ; then
    return
  fi
  
  aptcache=`which apt-cache`
  if [ -z $aptcache ] ; then
    return
  fi
  echo " ============= apt-get installation of ITPP =============== "
#  echo "apt-get and apt-cache were found!"
  verstring=`$aptcache -f search libitpp-dev | grep "Version"`
  aptitppver=(`echo "$verstring" | sed -e "s/Version\: //" -e "s/\./ /g"`)
  echo "libitpp-dev $verstring was in apt-cache detected"
  if [[ ${aptitppver[0]} -ne 4 ]]; then
    echo "But we need version 4... "
    echo
    return
  fi
  if [[ ${aptitppver[1]} -lt 2 ]]; then
    echo "But we need version 4.2... "
    echo
    return
  fi

  echo
  echo
  run_as_sudo "apt-get install libitpp-dev"
  
  dpkgstat=`dpkg -s libitpp-dev | grep -E "Status.* installed"`
  if [ ! -z "$dpkgstat" ]; then
    echo "Looks like itpp was installed successfully"
    itppinstalled=1
  fi
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

function install_blas {
  echo
  echo "BLAS not found. Trying to install BLAS..."
  echo
  echo "Press Enter to proceed,"
  echo "or Ctrl-C to break and install BLAS yourself."
  read

  download_file_with_forward http://www.netlib.org/blas/blas.tgz blas.tgz
  set -e
  tar -xzvf blas.tgz
  cd BLAS
  set +e
  if make ; then
    echo
    echo "BLAS Compiled! Installing..."
    echo
    run_as_sudo "cp blas_LINUX.a /usr/lib/libblas.a"
  else
    echo "BLAS Failed to compile..."
    echo "ITPP installation script cannot continue"
    exit 1
  fi
  cd ..
}

function install_lapack {
  echo
  echo "LAPACK not found. Trying to install LAPACK..."
  echo
  echo "Press Enter to proceed,"
  echo "or Ctrl-C to break and install LAPACK yourself."
  read
  download_file_with_forward http://www.netlib.org/lapack/lapack-3.3.1.tgz lapack.tgz
  set -e
  tar -xzvf lapack.tgz
  cd lapack-3.3.1
  cp make.inc.example make.inc
  echo -e "\nOPTS     = -O3" >> make.inc
  set +e
  if make lapacklib ; then
    echo
    echo "Lapack Compiled! Installing..."
    echo
    run_as_sudo "cp lapack_LINUX.a /usr/lib/liblapack.a"
  else
    echo "Lapack Failed to compile..."
    echo "ITPP installation script cannot continue"
    exit 1
  fi
  cd ..
}

function print_source_install_failure {
  echo "---------------------------------------------------"
  echo "| The ITPP Configuration cannot find BLAS/Lapack  |"
  echo "| and the gfortran compiler is not available.     |"
  echo "|                                                 |"
  echo "| This script will not be able to help you obtain |"
  echo "| all the required dependencies for ITPP          |"
  echo "|                                                 |"
  echo "| We suggest that you try to obtain binary        |"
  echo "| releases of BLAS/Lapack that are compatible     |"
  echo "| with your system and try again.                 |"
  echo "---------------------------------------------------" 
}

function install_from_source {
  if  [ $itppinstalled -eq 1 ]; then
    return
  fi
  
  echo " ============= Source installation of ITPP =============== "
  # detect wget
  echo "The script will now proceed to download itpp and try to install it"
  echo
  echo "Press Enter to proceed,"
  echo "or Ctrl-C to break and install itpp yourself, or try again with ITPP_ROOT"
  read
  echo "Downloading ITPP 4.2 from SourceForge..."
  blaslapack_install=$1
  download_file_with_forward http://sourceforge.net/projects/itpp/files/itpp/4.2.0/itpp-4.2.tar.gz/download itpp-4.2.tar.gz
  # unpack
  set -e
  tar -xzvf itpp-4.2.tar.gz
  cd itpp-4.2
  set +e
  blasinstalled=0
  lapackinstalled=0

  while true; do
    pushd . 
    f=`./configure --without-fft 2>&1`
    blascheck=`echo $f | grep "cannot find any BLAS library"`
    lapackcheck=`echo $f | grep "cannot find any LAPACK library"`
    if [[ ! -z $blascheck ]] || [[ ! -z $lapackcheck ]]; then
      if [[ -z `which gfortran` ]]; then
        print_source_install_failure
        exit 1
      else
        cd ..
        if [[ ! -z $blascheck ]]; then
          # if we have already tried to install blas
          if [ $blasinstalled -eq 1 ]; then
            print_source_install_failure
            exit 1
          fi
          install_blas
          blasinstalled=1
        fi
        if [[ ! -z $lapackcheck ]]; then
          if [ $lapackinstalled -eq 1 ]; then
            print_source_install_failure
            exit 1
          fi
          install_lapack
          lapackinstalled=1
        fi
      fi
    else
      popd
      break
    fi
    popd
  done
  # configure went well
  echo
  echo "ITPP Configuration Successful..."
  echo
  if make ; then
    echo
    echo "ITPP Compiled! Installing..."
    echo
    run_as_sudo "make install"
  else
    echo "ITPP Failed to compile..."
    echo "ITPP installation script cannot continue"
    exit 1
  fi
}

pushd .
mkdir -p $DEPS_DIR
cd $DEPS_DIR
echo "Detecting itpp..."
echo
# first argument 
test_for_itpp

if [ -z $itppfound ] ; then
  echo " ==================== ITPP Not Found ! ===================== "
  OS=`uname`
  itppinstalled=0
  if [[ "$OS" == "Darwin" ]]; then
    install_from_source
  elif [[ "$OS" == "Linux" ]]; then
    try_install_from_apt_get
    install_from_source
  else
    install_from_source
  fi
fi
find_itpp_lib
popd
