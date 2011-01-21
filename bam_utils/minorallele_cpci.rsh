#!/bin/bash
exec R --vanilla -q --slave -e "source(file=pipe(\"tail -n +4 $0\"))" --args $@
# this spawns R, reading this file from line #4...
# see: http://rwiki.sciviews.org/doku.php?id=tips:scriptingr

source(file="minorallele_cpci.R")
argv = commandArgs(trailingOnly=TRUE)
if (length(argv) == 3) {
    cat(CP.CI(as.integer(argv[1]),as.integer(argv[2]),as.integer(argv[3])));
    cat("\n");
} else {
    cat("Usage: minorallele_cpci.rsh n x num_alleles\n");
}
