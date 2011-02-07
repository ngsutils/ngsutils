#!/bin/bash

usage(){
    echo "Fetches dependencies required for working with BAM files"
    echo ""
    echo "Dependencies:"
    echo "  Cython"
    echo "  pysam"
    echo ""
    echo "Usage: `basename $0`"
    echo ""
    exit 1
}

if [ "$1" == "-h" ]; then
    usage
fi

if [ -e "`dirname $0`/ext/pysam" ]; then
    exit
fi

echo "Downloading and building dependencies..." >&2

WORK=`dirname $0`/ext/work

mkdir -p $WORK
cd $WORK
rm -rf *

python -c "import Cython" > /dev/null
if [ $? -ne 0 ]; then
    curl -LO "http://pypi.python.org/packages/source/C/Cython/Cython-0.14.1.tar.gz"
    tar zxvf Cython-0.14.1.tar.gz
    cd Cython-0.14.1
    python setup.py build
    python setup.py install --user
    cd ..
fi

curl -LO "http://pysam.googlecode.com/files/pysam-0.3.1.tar.gz"
tar zxvf pysam-0.3.1.tar.gz
cd pysam-0.3.1
python setup.py build

cp -r build/lib*/pysam ../..
cp build/lib*/*.so ../../pysam
