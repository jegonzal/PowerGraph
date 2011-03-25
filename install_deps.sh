#!/bin/bash
echo "===== install deps"
echo $1

source configure.deps


# see if BOOST_ROOT is the deps folder
# if it is, we may need to install it
if [ ! -z $BOOST_ROOT ] && [ $BOOST_ROOT == $PWD/deps ]; then
  echo
  echo "================================================================"
  echo "This script will now install Boost into"
  echo "$1/include and $1/lib"
  echo "However, only a small fraction of the Boost libraries were built"
  echo "to support the GraphLab library."
  echo "You may want not want to do this if you already have a complete"
  echo "installation of Boost elsewhere"
  echo
  echo "Press Y to continue, or N to skip the Boost installation"
  
  skipinstall=1
  while [ 1 ]; do
    read ans
    if [[ $ans == "N" || $ans == "n" ]]; then
      break
    elif [[ $ans == "Y" || $ans == "y" ]]; then
      skipinstall=0
      break
    else
      echo "Invalid Input: Press Y to continue, or N to skip the Boost installation"
    fi
   done
   if [ $skipinstall -eq 0 ]; then
    echo "Continuing Boost install"
    source ./boostinstall.sh
   else
    echo "Skipping Boost install"
   fi
fi