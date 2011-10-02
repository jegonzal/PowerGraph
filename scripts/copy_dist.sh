#!/bin/sh
cp dist/graphlabapi_v1_$1.tar.gz ~/www-graphlab/graphlab/release/
svn add ~/www-graphlab/graphlab/release/graphlabapi_v1_$1.tar.gz
svn commit -m "new release graphlabapi_v1_$1.tar.gz" ~/www-graphlab/graphlab
