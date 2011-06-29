#!/bin/sh

export CLASSPATH=$CLASSPATH:../build/:../lib/jython.jar

java -Xmx2000m  -Djava.library.path=../native/ demo.pagerank.PageRank ../python/testdata/stoch30x.csv

