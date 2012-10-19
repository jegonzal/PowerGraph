if [ ! -d src ]; then
  echo "Run from the graphlab root folder"
  exit
fi

pushd .
cd release/toolkits
make clean
popd

./configure --no_jvm -D NO_MPI:BOOL=true -D COMPILER_FLAGS="-mmacosx-version-min=10.5" -D MARCH=x86-64 -D MTUNE=generic
scripts/compile_static_release.sh

# now package a binary release
# for whatever reason the mac binaries are quite small... 
# stripping not necessary
rootdirname="graphlab_mac"
rm -rf ./$rootdirname
mkdir $rootdirname
for file in `cat scripts/binary_list.txt`
do
  dname=`dirname $file`
  mkdir -p $rootdirname/$dname
  cp release/$file $rootdirname/$dname/
done

#package all the rest of the stuff
#copy the license
mkdir $rootdirname/license
cp license/LICENSE.txt $rootdirname/license/

#copy the README
cp BINARY_README $rootdirname/
