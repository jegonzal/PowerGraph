#!/bin/sh
rm -f $1.reverse
head -n 1 $1 > $1.reverse
head -n 3 $1 | grep -v "%" | head -n 1 | awk '{print $1 " " $2 " " $3 }' >> $1.reverse
cat $1 | grep -v "%" | sort -n -k 2,2 -k 1,1 -T . >> $1
