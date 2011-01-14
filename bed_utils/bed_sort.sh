#!/bin/bash
#
# bed_sort.sh - Sorts BED files by chrom and start pos
#
# 2011-01-13 Marcus Breese
#
if [ "$1" == "" ]; then
    echo "Usage: `basename $0` in.bed [out.bed]"
    echo ""
    echo "Sorts a BED file.  If out.bed is given, the sorted output is"
    echo "written there.  Otherwise, in.bed is replaced with a sorted"
    echo "copy."
    echo ""
    exit 1
fi

if [ "$2" == "" ]; then
    TMP=".$1.sort.$$"
    sort -k1,1 -k2,2n $1 > $TMP
    mv $TMP $1
else
    sort -k1,1 -k2,2n $1 > $2
fi
