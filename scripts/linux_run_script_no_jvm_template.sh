#!/bin/bash
# This script will run the program in the same location and with the same
# name as this script (without the .sh). Passing it the same set of command
# line options.

# It also parses the following environment variables
# JAVA_HOME
#    Either JAVA_HOME or JVM_SO_PATH must be set.
#    This must point to the Java home directory.
#    For instance: /usr/lib/jvm/java-6-openjdk
#    This was tested with Oracle's implementation of Java (sun-jdk or open-jdk).

# JVM_SO_PATH
#    Either JAVA_HOME or JVM_SO_PATH must be set.
#    The script will expect to find libjvm.so in 
#    $JAVA_HOME/jre/lib/amd64/client/libjvm.so or
#    $JAVA_HOME/jre/lib/amd64/server/libjvm.so
#    If libjvm.so is not in either locations, the script will fail. In which 
#    case, you should set this variable to the directory containing libjvm.so.

# USE_SYSTEM_LIBS
#    Optional. If set to 1, the system's glibc (and other system dependencies) 
#   will be used instead of the provided versions. 




# get the program name. By convention we will make it so that the 
# script's name is the same as the program name. But with a ".sh" at the end
PROG=$0
#strip the ".sh" at the end of the script name
PROG=${PROG%.sh}

PROGDIR=`dirname $0`

JVM_PATH=""
HAS_HADOOP=0

if [ -z $USE_SYSTEM_LIBS ]; then
  USE_SYSTEM_LIBS=0
fi

if [ "$USE_SYSTEM_LIBS" -eq "1" ]; then
  echo "Using system libs."
  #using system libs
  $PROG $*
else
  # now. where do I find the dependency directory?
  # lets for now... assume that the directory organization must be
  # /gldeps
  # /toolkits/blah
  # /toolkits/otherblah 
  DEPDIR="$PROGDIR/../../gldeps"

  LIBPATH=$DEPDIR
  if [ ! -z "$JVM_PATH" ]; then
    LIBPATH=$LIBPATH:$JVM_PATH
  fi

  if [ ! -z "$LD_LIBRARY_PATH" ]; then
    LIBPATH=$LIBPATH:$LD_LIBRARY_PATH
  fi

  $DEPDIR/ld-linux-x86-64.so.2 --library-path $LIBPATH $PROG $*
fi

