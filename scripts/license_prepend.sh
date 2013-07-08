DIR="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $DIR
if ! grep -q Apache $1
then
echo $1
cat $DIR/../license/LICENSE_prepend.txt $1 > /tmp/out
mv /tmp/out $1
fi
