#!/bin/bash
REAL=`python -c 'import os,sys;print os.path.realpath(sys.argv[1])' "$0"`
DIR=`dirname "$REAL"`
SUBDIR=$(basename $0 | sed -e 's/utils//')

. "$DIR"/env/bin/activate
export PYTHONPATH=$PYTHONPATH:"$DIR"

if [ "$1" == "" ]; then
    # run all tests...
   	cd $DIR
    python -m unittest discover ngsutils 'test*py' -v
else
    coverage run $1
    if [ $? -eq 0 ]; then
        coverage report
    fi
fi
