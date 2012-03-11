#!/bin/sh
#script for converting matrix market format to its tranpose matrix
#written by danny bickson CMU
rm -f $1.reverse
head -n 1 $1 > $1.reverse
head -n 3 $1 | grep -v "%" | head -n 1 | awk '{print $2 " " $1 " " $3 }' >> $1.reverse
cat $1 | grep -v "%" | awk '{print $2 " " $1 " " $3}'| sort -n -k 2,2 -k 1,1 -T . >> $1.reverse
