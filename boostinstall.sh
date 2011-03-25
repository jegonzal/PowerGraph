pushd .
cd deps
if [ ! -d "boost_1_46_1" ] ; then
  echo "The boost installation in deps/ is missing! unable to continue!"
else
  cd boost_1_46_1
  # install with the new prefix
  echo "Now installing Boost... This could take a little while..."
  ./bjam --threading=multi --link=static --variant=release --prefix=$1 install
fi
popd