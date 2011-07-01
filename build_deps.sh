#!/bin/bash
CYTHON_VER=0.14.1
PYSAM_VER=0.4.1

usage(){
    echo "Fetches dependencies required for working with BAM files"
    echo ""
    echo "Dependencies:"
    echo "  Cython   (0.14.1)"
    echo "  pysam    (0.4.1)"
    echo ""
    echo "Usage: `basename $0`"
    echo ""
    exit 1
}

control_c() {
    rm -rf "`dirname $0`/work/*"
    exit 1
}

if [ "$1" == "-h" ]; then
    usage
fi

if [ "`which curl`" == "" ]; then
    echo "Missing dependency: curl"
    exit 1
fi

PYMAJOR=`python -V 2>&1 | awk '{print $2}' | sed -e 's/\./ /g' | awk '{print $1}'`
PYMINOR=`python -V 2>&1 | awk '{print $2}' | sed -e 's/\./ /g' | awk '{print $2}'`

if [ "$PYMAJOR" != "2" ]; then
    echo "ngsutils requires Python version 2.6+. It has not been tested with Python 3."
    exit 0
fi
if [ $PYMINOR -lt 6 ]; then
    echo "ngsutils requires Python version 2.6+. It has not been tested with Python 3."
    exit 0
fi

if [ -e "`dirname $0`/ext/pysam" ]; then
    exit 0
fi

trap control_c SIGINT

echo "Downloading and building dependencies..." >&2

ABSPATH=`python -c 'import os,sys;print os.path.realpath(sys.argv[1])' "$0"`
WORK=`dirname $ABSPATH`/ext/work
mkdir -p "$WORK"
cd "$WORK"
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
curl -LO "http://pysam.googlecode.com/files/pysam-0.4.1.tar.gz"
tar zxf pysam-0.4.1.tar.gz
cd pysam-0.4.1

echo "  [pysam] building" >&2
python setup.py build 2>> $WORK/build.log >> $WORK/build.log
if [ $? -ne 0 ]; then
    echo "  [pysam] error building - see work/build.log" >&2
    exit 1
fi

cp -r build/lib*/pysam ../..
cp build/lib*/*.so ../../pysam
echo "  Done!" >&2
echo "" >&2

