#!/usr/bin/env python
## category General
## desc Extract genomic/transcript sizes for genes
'''
Extracts the genomic and transcript sizes for genes in GTF format
'''

import sys
import os
from ngsutils.gtf import GTF


def usage():
    print __doc__
    print '''\
Usage: gtfutils genesize filename.gtf
'''
    sys.exit(1)


def gtf_genesize(gtf, out=sys.stdout):
    out.write('#gene\tgenomic-size\ttranscript-size\n')
    for gene in gtf.genes:
        cols = [gene.gene_name]
        cols.append((gene.end - gene.start))

        maxsize = 0
        for txs in gene.transcripts:
            size = 0
            for start, end in txs.exons:
                size += (end - start)
            maxsize = max(size, maxsize)

        cols.append(maxsize)
        out.write('%s\n' % '\t'.join([str(x) for x in cols]))


if __name__ == '__main__':
    fname = None
    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if not fname and os.path.exists(arg):
            fname = arg
        else:
            usage()

    if not fname:
        usage()

    gtf = GTF(fname)
    gtf_genesize(gtf)
