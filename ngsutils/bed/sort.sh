#!/bin/bash
## category General
## desc Sorts a BED file (in place)
#
# bed_sort.sh - Sorts BED files by chrom and start pos
#
# 2011-01-13 Marcus Breese
#
if [[ "$1" == "" || "$1" == "-h" ]]; then
    echo "Usage: bedutils sort in.bed [out.bed]"
    echo ""
    echo "Sorts a BED file.  If out.bed is given, the sorted output is"
    echo "written there.  Otherwise, in.bed is replaced with a sorted"
    echo "copy."
    echo ""
    exit 1
fi

if [ "$1" == "-" ]; then
    sort -k1,1 -k2,3n -
elif [ "$2" == "" ]; then
    DN=`dirname "$1"`
    BASE=`basename "$1"`
    TMP="$DN/.$BASE.sort.$$"
    sort -k1,1 -k2,3n "$1" > "$TMP"
    mv "$TMP" "$1"
else
    sort -k1,1 -k2,3n "$1" > "$2"
fi
