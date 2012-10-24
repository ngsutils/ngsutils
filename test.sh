#!/bin/bash

if [ "$1" == "" ]; then
    # run all tests...
    python -m unittest discover ngsutils 'test*py' -v
else
    PYTHONPATH=. coverage run $1
    if [ $? -eq 0 ]; then
        coverage report
    fi
fi
