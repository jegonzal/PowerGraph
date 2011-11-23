#!/bin/sh
cp dist/graphlabapi_v1_$1.tar.gz ~/www-graphlab/graphlab/release/
svn add ~/www-graphlab/graphlab/release/graphlabapi_v1_$1.tar.gz
vim ~/www-graphlab/graphlab/download.html
cd ~/www-graphlab/graphlab/
svn commit -m "new release graphlabapi_v1_$1.tar.gz" 

