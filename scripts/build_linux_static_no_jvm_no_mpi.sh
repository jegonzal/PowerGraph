if [ ! -d src ]; then
  echo "Run from the graphlab root folder"
  exit
fi
./configure -D MARCH=x86-64 -D MTUNE=generic --no_jvm -D NO_MPI:BOOL=true -D COMPILER_FLAGS:STRING="-static-libgcc\ -static-libstdc++" 
scripts/compile_static_release.sh $@

# is this a openmpi or a mpich2 build?
rootdirname="graphlab_no_jvm_no_mpi"
unstrippeddirname="graphlab_unstripped_no_jvm_no_mpi"
ISOPENMPI=0


# now package a binary release
rm -rf ./$rootdirname
rm -rf ./$unstrippeddirname
mkdir $rootdirname
mkdir $unstrippeddirname
mkdir $rootdirname/gldeps
mkdir $unstrippeddirname/gldeps

tmp=$@
if test $# -lt 1 ; then
  tmp=`cat scripts/binary_list.txt`
fi


for file in $tmp
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

#pack
tar -cjvf $rootdirname.tar.bz2 $rootdirname
tar -cjvf $unstrippeddirname.tar.bz2 $unstrippeddirname
