pushd .
cd deps
if [ ! -d "kyotocabinet-1.2.53" ] ; then
  echo "The Kyoto Cabinet installation in deps/ is missing! unable to continue!"
else
  cd kyotocabinet-1.2.53
  # install with the new prefix
  echo "Now installing Kyoto Cabinet... "
  ./configure --prefix=$installprefix 
  make install
fi
popd
