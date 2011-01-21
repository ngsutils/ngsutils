#!/bin/bash
exec R --vanilla -q --slave -e "source(file=pipe(\"tail -n +4 $0\"))" --args $@
# this spawns R, reading this file from line #4...
# see: http://rwiki.sciviews.org/doku.php?id=tips:scriptingr

source(file="minorallele_cpci.R")
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

argv = commandArgs(trailingOnly=TRUE)
if (length(argv) == 1) {
    close_socket <- FALSE;
    conn <- make.socket(host = "localhost", as.integer(argv[1]), fail = TRUE, server = TRUE);
    on.exit(close.socket(conn));
    
    while (!close_socket) {
        line <- read.socket(conn,1024,TRUE);
        cols <- strsplit(trim(line),' ');
        if (cols[[1]][1] == "quit") {
            close_socket = TRUE;
        } else if (length(cols[[1]]) == 3) {
            val <- CP.CI(as.integer(cols[[1]][1]),as.integer(cols[[1]][2]),as.integer(cols[[1]][3]));
            write.socket(conn,toString(val[1]));
            write.socket(conn," ");
            write.socket(conn,toString(val[2]));
            write.socket(conn,"\n");
        }
    }
    
} else if (length(argv) == 3) {
    cat(CP.CI(as.integer(argv[1]),as.integer(argv[2]),as.integer(argv[3])));
    cat("\n");
} else {
    cat("Usage: minorallele_cpci.rsh n x num_alleles\n");
    cat("       minorallele_cpci.rsh port (sets up a port for repeated calculations)\n");
}
