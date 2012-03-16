#!/bin/sh
rm -f $1.out
head -n 1 $1 > $1.out
head -n 3 $1 | grep -v "%" | head -n 1 | awk '{print $1 " " $2 " " $3 }' >> $1.out
cat $1 | grep -v "%" | sort -n -k 1,1 -k 2,2 -T . >> $1.out
