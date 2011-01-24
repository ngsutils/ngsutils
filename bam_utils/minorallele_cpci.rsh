#!/bin/bash
exec R --vanilla -q --slave -e "source(file=\"$(dirname $0)/minorallele_cpci.R\");source(file=pipe(\"tail -n +4 $0\"))" --args $@
# this spawns R, reading this file from line #4...
# see: http://rwiki.sciviews.org/doku.php?id=tips:scriptingr

argv = commandArgs(trailingOnly=TRUE)
if (length(argv) == 2) {
    alleles <- as.integer(argv[1]);
    con <- file(argv[2]);
    open(con);
    
    #header
    line <- readLines(con, 1);
    cat(line);
    
    cat("\tMean level\t95% CI low\t95% CI high\tCI Range\tlow count\thigh count\n");
    
    while (length(line <- readLines(con,1,warn=FALSE)) > 0) {
        cols <- strsplit(line,'\t')[[1]];

        ref <- as.integer(cols[9])
        alt <- as.integer(cols[10])
        
        val <- CP.CI(ref+alt,alt,alleles);
        ci_high = val[2];
        ci_low = val[1];

        cols[11] <- alt / (alt+ref);
        cols[12] <- ci_low;
        cols[13] <- ci_high;
        cols[14] <- ci_high-ci_low;
        cols[15] <- ci_low * alleles;
        cols[16] <- ci_high * alleles;
        
        i <- 1;
        while (i<17) {
            if (i > 1) {
                cat("\t");
            }
            cat(cols[i]);
            i <- i + 1
        }
        cat("\n");
    }
    close(con);
    
} else if (length(argv) == 3) {
    cat(CP.CI(as.integer(argv[1]),as.integer(argv[2]),as.integer(argv[3])));
    cat("\n");
} else {
    cat("Usage: minorallele_cpci.rsh n x num_alleles\n");
    cat("       minorallele_cpci.rsh num_alleles input.txt > output\n");
}
