if [ ! -d src ]; then
  echo "Run from the graphlab root folder"
  exit
fi
./configure -D MARCH=x86-64 -D MTUNE=generic -D COMPILER_FLAGS="-static-libgcc\ -static-libstdc++"
scripts/compile_static_release.sh

# is this a openmpi or a mpich2 build?
ISOPENMPI=0
if grep -q mpi_cxx config.log
then
  rootdirname="graphlab_openmpi"
  unstrippeddirname="graphlab_openmpi_unstripped"
  ISOPENMPI=1
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

  deps=$(ldd release/$file | awk 'BEGIN{ORS=" "}$1 \
      ~/^\//{print $1}$3~/^\//{print $3}' \
       | sed 's/,$/\n/')

  for dep in $deps
  do
    depname=`basename $dep`
    # definitely exclude jvm
    if [[ $depname == "libjvm.so" ]]; then
      continue
    fi
    if [ ! -f "$rootdirname/gldeps/$depname" ]; then
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
cp BINARY_README $rootdirname/README
cp BINARY_README $unstrippeddirname/README

# I am unable to get openmpi to work properly with the ld hack
# since it appears to have complicated binary dependencies. 
# (it forks and launches some other daemon which has its own dependencies)
# I will give up on this for now and try to get ABI compatibility.
# it seems like 1.3 is compatible with 1.4 and 1.5 is compatbile with 1.6
if [ ISOPENMPI -eq 1 ]; then
  rm $rootdirname/gldeps/libmpi.* $rootdirname/gldeps/libopen-*
  rm $unstrippeddirname/gldeps/libmpi.* $unstrippeddirname/gldeps/libopen-*
fi

#pack
tar -cjvf $rootdirname.tar.bz2 $rootdirname
tar -cjvf $unstrippeddirname.tar.bz2 $unstrippeddirname
