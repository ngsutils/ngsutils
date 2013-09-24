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
    out.write('#gene\tgenomic-size\ttranscript-size\ttotal-intron-length\tcoding-length\t5\'UTR-length\t3\'UTR-length\n')
    for gene in gtf.genes:
        cols = [gene.gene_name]
        cols.append((gene.end - gene.start))

        maxtx = 0
        maxintron = 0
        maxcoding = 0
        maxutr5 = 0
        maxutr3 = 0

        for txs in gene.transcripts:
            size = 0
            intron_size = 0
            last_end = None
            for start, end in txs.exons:
                size += (end - start)

                if last_end is not None:
                    intron_size += (start - last_end)
                last_end = end

            cds = 0
            utr5 = 0
            utr3 = 0

            if txs.has_cds:
                for s, e in txs.cds:
                    cds += (e - s)

                for s, e in txs.utr_5:
                    utr5 += (e - s)

                for s, e in txs.utr_3:
                    utr3 += (e - s)



            maxtx = max(size, maxtx)
            maxintron = max(intron_size, maxintron)
            maxcoding = max(cds, maxcoding)
            maxutr5 = max(utr5, maxutr5)
            maxutr3 = max(utr3, maxutr3)

        cols.append(maxtx)
        cols.append(maxintron)
        cols.append(maxcoding)
        cols.append(maxutr5)
        cols.append(maxutr3)
        
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
