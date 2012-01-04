#!/bin/sh

export CLASSPATH=$CLASSPATH:../build/:../lib/jython.jar

cd ../python
java -Xmx2000m -Djava.library.path=../native/ graphlab.wrapper.PythonGraphlab lasso_load.py lasso_config.py lasso_update.py lasso_post.py testdata/testphoto.txt

