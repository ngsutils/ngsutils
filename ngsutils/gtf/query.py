#!/usr/bin/env python
## category General
## desc Query a GTF file by coordinates
'''
Query a GTF file by coordinates

This will return the gene (or genes) present in a given genome region.
'''

import sys
import os
from ngsutils.gtf import GTF


def usage(msg=None):
    if msg:
        print '%s\n' % msg
    print __doc__
    print 'Usage: gtfutils query filename.gtf{.gz} chrom:start-end'
    sys.exit(1)


def gtf_query(fname, chrom, start, end, strand=None):
    gtf = GTF(fname)
    for gene in gtf.find(chrom, start, end, strand):
        yield gene

if __name__ == '__main__':
    filename = None
    region = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif not filename and os.path.exists(arg):
            filename = arg
        elif not region:
            region = arg

    if not filename:
        usage('Missing input file')
    if not region:
        usage('Missing query region! (chrom:start-end)')

    chrom, startend = region.split(':')
    start, end = [int(x) for x in startend.split('-', 1)]

    for gene in gtf_query(filename, chrom, start, end):
        print  ' %s %s:%s-%s[%s]' % (gene.gene_name, gene.chrom, gene.start, gene.end, gene.strand)
