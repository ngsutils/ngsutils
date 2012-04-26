#!/bin/bash

if [ ! -e env ]; then 
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

    echo "Initializing virtualenv folder (env)"
    $VIRTUALENV env
fi

. env/bin/activate

if [ "$1" == "-q" ]; then
    pip install -qr requirements.txt
else
echo "Installing required libraries"
    pip install -r requirements.txt
fi
