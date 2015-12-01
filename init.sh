#!/bin/sh

if [ "$PYTHON" = "" ]; then
    PYTHON="python"
fi

# Use embedded virtualenv
REAL=`python -c 'import os,sys;print os.path.realpath(sys.argv[1])' "$0"`
VIRTUALENV="$PYTHON $(dirname $REAL)/support/virtualenv.py"

if [ ! -e venv ]; then 
    echo "Initializing virtualenv folder (venv)"
    $VIRTUALENV --no-site-packages venv
fi

. venv/bin/activate

if [ $(uname -s) = "Darwin" ]; then
    # Mac OS X Mountain Lion compiles with clang by default...
    # clang and cython don't get along... so force it to use gcc

    if [ "$(cc --version | grep -i clang)" != "" ]; then
        echo "Using GCC"
        export CC=/usr/bin/gcc
        export CXX=/usr/bin/g++
    fi
fi

PYTHONMAJOR=$($PYTHON -V 2>&1 | sed -e 's/\./ /g' | awk '{print $2}')
PYTHONMINOR=$($PYTHON -V 2>&1 | sed -e 's/\./ /g' | awk '{print $3}')

if [ "$PYTHONMAJOR" -ne 2 ]; then
    echo "Requires Python 2.6+"
    exit
fi
if [ "$PYTHONMINOR" -lt 6 ]; then
    echo "Requires Python 2.6+"
    exit
fi
if [ "$PYTHONMINOR" -eq 6 ]; then
    pip install unittest2
fi

echo "Installing required libraries"
#pip install cython==0.16
pip install -r requirements.txt
