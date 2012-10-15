if [ ! -d src ]; then
  echo "Run from the graphlab root folder"
else
  ./configure --no_jvm -D NO_MPI:BOOL=true -D COMPILER_FLAGS="-mmacosx-version-min=10.4" -D MARCH=x86-64 -D MTUNE=generic
  scripts/build_static.sh
fi
