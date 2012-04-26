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

echo "Installing required libraries"
pip install -r requirements.txt
