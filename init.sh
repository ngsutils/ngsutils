#!/bin/bash

VIRTUALENV=`which virtualenv-2.7`
if [ "$VIRTUALENV" == "" ]; then
    VIRTUALENV=`which virtualenv-2.6`
fi
if [ "$VIRTUALENV" == "" ]; then
    VIRTUALENV=`which virtualenv`
fi
if [ "$VIRTUALENV" == "" ]; then
    echo "Missing virtualenv!"
    exit 1
fi

if [ ! -e env ]; then 
    echo "Initializing virtualenv folder (env)"
    $VIRTUALENV --no-site-packages env
fi

. env/bin/activate

if [ $(uname -s) == "Darwin" ]; then
    # Mac OS X Mountain Lion compiles with clang by default...
    # clang and cython don't get along... so force it to use gcc

    if [ "$(cc --version | grep -i clang)" != "" ]; then
        echo "Using GCC"
        export CC=/usr/bin/gcc
        export CXX=/usr/bin/g++
    fi
fi

echo "Installing required libraries"
pip install -r requirements.txt
