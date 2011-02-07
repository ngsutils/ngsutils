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

control_c() {
    rm -rf `dirname $0`/ext/*
    exit 1
}

if [ "$1" == "-h" ]; then
    usage
fi

if [ -e "`dirname $0`/ext/pysam" ]; then
    exit
fi

trap control_c SIGINT

echo "Downloading and building dependencies..." >&2

WORK=`dirname $0`/ext/work
mkdir -p $WORK
cd $WORK
touch build.log
rm -rf *

python -c "import Cython" &> /dev/null
if [ $? -ne 0 ]; then
    echo "  [cython] downloading" >&2
    curl -LO "http://pypi.python.org/packages/source/C/Cython/Cython-0.14.1.tar.gz"
    tar zxf Cython-0.14.1.tar.gz
    cd Cython-0.14.1

    echo "  [cython] building" >&2
    python setup.py build 2>> $WORK/build.log >> $WORK/build.log
    if [ $? -ne 0 ]; then
        echo "  [cython] error building - see work/build.log" >&2
        exit 1
    fi

    echo "  [cython] installing" >&2
    python setup.py install --user 2>> $WORK/build.log >> $WORK/build.log
    if [ $? -ne 0 ]; then
        echo "  [cython] error installing - see work/build.log" >&2
        exit 1
    fi
    cd ..
fi

echo "  [pysam] downloading" >&2
curl -LO "http://pysam.googlecode.com/files/pysam-0.3.1.tar.gz"
tar zxf pysam-0.3.1.tar.gz
cd pysam-0.3.1

echo "  [pysam] building" >&2
alias gcc='gcc -D_GNU_SOURCE'
python setup.py build 2>> $WORK/build.log >> $WORK/build.log
unalias gcc
if [ $? -ne 0 ]; then
    echo "  [pysam] error building - see work/build.log" >&2
    exit 1
fi

cp -r build/lib*/pysam ../..
cp build/lib*/*.so ../../pysam
echo "  Done!" >&2
echo "" >&2

