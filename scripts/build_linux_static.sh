if [ ! -d src ]; then
  echo "Run from the graphlab root folder"
  exit
fi
pushd .
cd release/toolkits
make clean
popd
./configure -D MARCH=x86-64 -D MTUNE=generic
scripts/compile_static_release.sh

# is this a openmpi or a mpich2 build?
if grep -q openmpi config.log
then
  rootdirname="graphlab_openmpi"
  unstrippeddirname="graphlab_openmpi_unstripped"
elif grep -q mpich config.log
then
  rootdirname="graphlab_mpich2"
  unstrippeddirname="graphlab_mpich2_unstripped"
else
  echo "Unable to detect MPI type"
  exit
fi


# now package a binary release
rm -rf ./$rootdirname
rm -rf ./$unstrippeddirname
mkdir $rootdirname
mkdir $unstrippeddirname
mkdir $rootdirname/gldeps
mkdir $unstrippeddirname/gldeps

for file in `cat scripts/binary_list.txt`
do
  dname=`dirname $file`
  fname=`basename $file`

  deps=$(ldd $file | awk 'BEGIN{ORS=" "}$1 \
      ~/^\//{print $1}$3~/^\//{print $3}' \
       | sed 's/,$/\n/')

  for dep in $deps
  do
    depname=`basename $dep`
    if [ ! -a "$rootdirname/gldeps/$depname" ]; then
      echo "Copying $dep"
      cp "$dep" "$rootdirname/gldeps/"
      cp "$dep" "$unstrippeddirname/gldeps/"
    fi
  done

  mkdir -p $rootdirname/$dname
  cp release/$file $rootdirname/$dname/
  #strip it
  strip $rootdirname/$dname/$fname
  #package the script
  cp scripts/linux_run_script_template.sh $rootdirname/$dname/$fname.sh

  #repeat for unstripped
  mkdir -p $unstrippeddirname/$dname
  cp release/$file $unstrippeddirname/$dname/
  #package the script
  cp scripts/linux_run_script_template.sh $unstrippeddirname/$dname/$fname.sh
done

#package all the rest of the stuff
#copy the license
mkdir $rootdirname/license
cp license/LICENSE.txt $rootdirname/license/

mkdir $unstrippeddirname/license
cp license/LICENSE.txt $unstrippeddirname/license/

#copy the README
cp BINARY_README $rootdirname/
cp BINARY_README $unstrippeddirname/
