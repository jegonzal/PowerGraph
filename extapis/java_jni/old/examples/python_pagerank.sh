#!/bin/sh

export CLASSPATH=$CLASSPATH:../build/:../lib/jython.jar

cd ../python
java -Xmx2000m -Djava.library.path=../native/ graphlab.wrapper.PythonGraphlab pagerank_load.py pagerank_config.py pagerank_update.py pagerank_post.py testdata/stoch30x.csv


