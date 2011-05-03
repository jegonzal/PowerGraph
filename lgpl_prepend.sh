tmpfile=`mktemp`
for f in `find src/graphlab/*  \( \( -name "*.cpp" -or -name "*.hpp" \) -not -wholename "src/graphlab/extern/*" \)`
do
  if ! grep -q "Lesser General Public License" $f 
  then
    echo $f
    cat LGPL_prepend.txt $f > $tmpfile
    mv $tmpfile $f
  fi
done