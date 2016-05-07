#!/bin/bash
REAL=`python -c 'import os,sys;print os.path.realpath(sys.argv[1])' "$0"`
DIR=`dirname "$REAL"`
SUBDIR=$(basename $0 | sed -e 's/utils//')

. "$DIR"/venv/bin/activate
export PYTHONPATH=$PYTHONPATH:"$DIR"
export HIDE_ETA="1"
export TESTING="1"

if [ "$1" = "" ]; then
    # run all tests...
    cd $DIR
    if [ "`which unit2`" != "" ]; then
        unit2 discover ngsutils 'test_*py' -v
    else
        python -m unittest2 discover ngsutils 'test_*py' -v
    fi
else
    while [ "$1" != "" ]; do
        echo "$1"
        coverage run $1
        if [ $? -eq 0 ]; then
            coverage report
        fi
        shift
    done
fi
