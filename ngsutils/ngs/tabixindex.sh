#!/bin/bash

usage() {
	echo "Tabix index a text file"
	echo ""
	echo "This script automates the process of compressing a text file with 'bgzip' and"
	echo "indexing the file using 'tabix'. Each of these programs must be pre-installed"
	echo "and part of your path."
	echo ""
	echo "Usage: ngsutils tabixindex file.txt{.gz} {tabix options}"
	echo ""
	echo "Tabix options:"
	START="-n+$(tabix 2>&1 | grep -n 'Options:' | sed -e 's/\:/ /' | awk '{print $1}')"
	tabix 2>&1 | tail $START
	exit 1
}

if [[ "$(which tabix)" == "" || "$(which bgzip)" == "" ]]; then
	echo "tabix and bgzip must be in your PATH"
	usage
fi

if [ "$1" == "" ]; then
	echo "Missing filename!"
	usage
fi

FNAME="$1"
shift

if [ "$1" == "" ]; then
	echo "Missing Tabix options!"
	usage
fi

DIR=$(dirname "$FNAME")
NEWNAME="$(basename "$FNAME" | sed -e 's/\.gz$//').bgz"

# use pv, if available...
CAT=cat
if [ "$(which pv)" != "" ]; then
	CAT=pv
fi

if [ ! -e "$DIR/$NEWNAME" ]; then
	if [ "$(basename "$FNAME" | sed -e 's/\.gz$//')" == "$(basename "$FNAME")" ]; then
		# not gzipped
		echo "Compressing $FNAME with bgzip"
		$CAT $FNAME | bgzip > "$DIR/.$NEWNAME.tmp"
		mv "$DIR/.$NEWNAME.tmp" "$DIR/$NEWNAME"
	else
		# gzipped source
		echo "Re-compressing $FNAME with bgzip"
		$CAT $FNAME | gunzip | bgzip > "$DIR/.$NEWNAME.tmp"
		mv "$DIR/.$NEWNAME.tmp" "$DIR/$NEWNAME"
	fi
fi

echo "Indexing..."
tabix "$@" $DIR/$NEWNAME
if [ $? -ne 0 ]; then
	rm $DIR/$NEWNAME
fi
