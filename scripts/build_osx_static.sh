if [ ! -d src ]; then
  echo "Run from the graphlab root folder"
  exit
fi


./configure --no_jvm -D NO_MPI:BOOL=true -D COMPILER_FLAGS="-mmacosx-version-min=10.7" -D MARCH=x86-64 -D MTUNE=generic -D HAS_CRC32:BOOL=false
scripts/compile_static_release.sh

echo "Packaging binary release..."

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
cp BINARY_README $rootdirname/README

echo "Binary release packaged in $rootdirname"
tar -cjvf $rootdirname.tar.bz2 $rootdirname
