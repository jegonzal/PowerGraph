#!/bin/bash

major_version=2.1

echo "THIS MUST BE RUN IN GRAPHLAB HOME"

## JOEY: WHY ARE WE REMOVING THE FOLDER AND THEN USING RSYNC?
rm -fR dist/graphlabapi
mkdir -p dist/graphlabapi
rsync -vv -al --delete --delete-excluded \
    --exclude=/debug --exclude=/release --exclude=/profile --exclude=/apps \
    --exclude=.hg --exclude=/matlab \
    --exclude=/dist --exclude=/deps --exclude=*~ --exclude=*.orig --exclude=/configure.deps \
    --exclude /make_dist --exclude /BINARY_README * dist/graphlabapi/.

mkdir dist/graphlabapi/apps
cp dist/graphlabapi/demoapps/CMakeLists.txt dist/graphlabapi/apps/
version=`hg summary | grep parent | sed 's/parent: //g' | sed 's/:.*//g'`
version="v${major_version}.$version"
echo "Version: $version"


cd dist
tar -vz \
    -cf graphlabapi_${version}.tar.gz \
    graphlabapi
cd ..

ls -al dist | tail -n 1

