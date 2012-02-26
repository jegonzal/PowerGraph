#!/bin/bash

pushd .
mkdir ./deps
cd ./deps


function test_for_hadoop {
  # if force hadoop install is set, don't probe. assume failure
    if [ ! -z $force_hadoop_install ] ; then
        return
    fi
    if [ ! -z $HADOOP_HOME ] ; then
        echo "Hadoop home found at:  $HADOOP_HOME"
        hadoopfound=1
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



echo "Probing for hadoop..."
test_for_hadoop
if [ ! -z $force_hadoop_install ] || [ -z $hadoopfound ]; then
    echo " ==================== Hadoop Not Found ! ==================== "
    echo
    echo "Hadoop libraries were not detected. Either your system does not have"
    echo "Hadoop, or HADOOP_HOME was not specified correctly."
    echo "The script will now proceed to download hadoop and install it locally"
    echo "in the graphlabapi/deps directory."
    echo
    echo "Press Enter to proceed,"
    echo "or Ctrl-C to break and install Hadoop yourself, or try again with HADOOP_HOME"
    if [ -z $always_yes ] ; then
        read
    fi
#    hadoop_ver=hadoop-1.0.0
    hadoop_ver=hadoop-0.20.2
    echo "Download $hadoop_ver"
    #download hadoop
    rm $hadoop_ver.tar.gz
    download_file_with_forward \
        http://mirror.cc.columbia.edu/pub/software/apache//hadoop/common/$hadoop_ver/$hadoop_ver.tar.gz \
        $hadoop_ver.tar.gz
    echo "Cleaning up prior installations of hadoop..."
    if [ -d $hadoop_ver ] ; then
        rm -rf $hadoop_ver
    fi
    tar -xzvf $hadoop_ver.tar.gz
    # set HADOOP_HOME
    HADOOP_HOME="$PWD/$hadoop_ver"
    echo "HADOOP_HOME=$HADOOP_HOME" >> ../configure.deps
else
    echo "Hadoop Found!"
fi

popd



