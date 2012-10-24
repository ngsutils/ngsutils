#!/bin/bash

if [ "$1" == "" ]; then
    # run all tests...
    python -m unittest discover ngsutils 'test*py' -v
else
    PYTHONPATH=. python $1 -v
fi
