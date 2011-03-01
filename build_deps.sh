#!/bin/bash

usage(){
    echo "Fetches dependencies required for working with BAM files"
    echo ""
    echo "Dependencies:"
    echo "  Cython   (0.14.1)"
    echo "  pysam    (0.3.1)"
    echo "  tabix    (tabix-0.2.3)"
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

if [ "`which curl`" == "" ]; then
    echo "Missing dependency: curl"
    exit 1
fi

if [ "`which tabix`" == "" ]; then
    mkdir -p ~/bin/src
    cd ~/bin/src

    if [ ! -e tabix-0.2.3 ]; then
        echo "[tabix] downloading"
        curl -LO "http://downloads.sourceforge.net/project/samtools/tabix/tabix-0.2.3.tar.bz2"
        tar jxvf tabix-0.2.3.tar.bz2
    fi

    cd tabix-0.2.3
    echo "[tabix] building"
    make
    cd ../..
    if [ -e tabix ]; then
        rm tabix
    fi
    if [ -e bgzip ]; then
        rm bgzip
    fi
    ln -s src/tabix-0.2.3/tabix .
    ln -s src/tabix-0.2.3/bgzip .
    cd src

    FOUND="0"
    for f in `echo $PATH | sed -e 's/:/ /g'`; do
        if [[ "$f" == "$HOME/bin" || "$f" == "~/bin" ]]; then
            FOUND="1"
        fi
    done

    if [ "$FOUND" == "0" ]; then
        echo "export PATH=$PATH:$HOME/bin" >> $HOME/.bashrc
        . $HOME/.bashrc
    fi
fi    

if [ -e "`dirname $0`/ext/pysam" ]; then
    exit 0
fi

trap control_c SIGINT

echo "Downloading and building dependencies..." >&2

ABSPATH=`python -c 'import os,sys;print os.path.realpath(sys.argv[1])' "$0"`
WORK=`dirname $ABSPATH`/ext/work
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
python setup.py build 2>> $WORK/build.log >> $WORK/build.log
if [ $? -ne 0 ]; then
    echo "  [pysam] error building - see work/build.log" >&2
    exit 1
fi

cp -r build/lib*/pysam ../..
cp build/lib*/*.so ../../pysam
echo "  Done!" >&2
echo "" >&2

