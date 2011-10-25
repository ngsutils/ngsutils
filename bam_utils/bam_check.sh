#!/bin/bash

if [ "$1" == "" ]; then
    echo "Usage: bamutils check filename.bam"
    exit -1
fi

if [ "`which samtools`" == "" ]; then
    echo 'Missing samtools from $PATH'
    exit -1
fi
echo -n "$1 "
samtools view "$1" &> /dev/null
if [ $? -eq 0 ]; then
    echo "OK"
    exit 0
else
    echo "ERROR"
    exit 1
fi
