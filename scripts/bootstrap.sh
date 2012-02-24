#!/bin/bash

IN_BOOTSTRAP=1
DEPS_DIR=deps
mkdir -p $DEPS_DIR

source scripts/cmakecheck.sh
echo "CMAKE=$cmakecmd" >> configure.deps


source scripts/boostcheck.sh
# we add a boost root only if it is the deps directory
if [ ! -z $BOOST_ROOT ] && [ $BOOST_ROOT == $PWD/$DEPS_DIR ]; then
  echo "BOOST_ROOT=$BOOST_ROOT" >> configure.deps
fi

source scripts/kccheck.sh
# we add a boost root only if it is the deps directory
if [ ! -z $KC_ROOT ] && [ $KC_ROOT == $PWD/$DEPS_DIR ]; then
  echo "KC_ROOT=$KC_ROOT" >> configure.deps
fi


# if [ ! -z $YRL_EXPERIMENTAL ]; then
#   echo "Installing AVRO"
#   source scripts/avrocheck.sh
#   if [ ! -z $AVRO_ROOT ] && [ $AVRO_ROOT == $PWD/$DEPS_DIR ]; then
#     echo "AVRO_ROOT=$AVRO_ROOT" >> configure.deps
#   fi
# fi

