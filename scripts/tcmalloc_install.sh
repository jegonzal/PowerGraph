#!/bin/bash

pushd .
mkdir ./deps
cd ./deps

install_location=$PWD


function test_for_tcmalloc {
  # if force tcmalloc install is set, don't probe. assume failure
    if [ ! -z $force_tcmalloc_install ] ; then
        return
    fi
    if [ ! -z $TCMALLOC_ROOT ] ; then
        echo "Tcmalloc home found at:  $TCMALLOC_ROOT"
        tcmallocfound=1
    fi
}



echo "Probing for tcmalloc..."
test_for_tcmalloc
if [ ! -z $force_tcmalloc_install ] || [ -z $tcmallocfound ]; then
    echo " ==================== Tcmalloc Not Found ! ==================== "
    echo
    echo "Tcmalloc libraries were not detected. Either your system does not have"
    echo "Tcmalloc, or TCMALLOC_ROOT was not specified correctly."
    echo "The script will now proceed to download tcmalloc and install it locally"
    echo "in the graphlabapi/deps directory."
    echo
    echo "Press Enter to proceed,"
    echo "or Ctrl-C to break and install Tcmalloc yourself, or try again with TCMALLOC_ROOT"
    if [ -z $always_yes ] ; then
        read
    fi
    echo "Download TCMalloc"
    svn checkout http://gperftools.googlecode.com/svn/trunk/ gperftools
    cd gperftools
    ./configure --enable-frame-pointers --prefix $install_location
    make -j2
    make install
    cd ..
    # set TCMALLOC_ROOT
    TCMALLOC_ROOT=$install_location
    echo "TCMALLOC_ROOT=$TCMALLOC_ROOT" >> ../configure.deps
else
    echo "Tcmalloc Found!"
fi

popd






