#!/usr/bin/env python
## category General
## desc Builds a TopHat-compatible .juncs file based on GFF/GTF model
'''
Takes a GTF gene model outputs a TopHat-compatible .juncs file.

This file is a tab-delimited file with the following fields:
chrom    3'end     5'start    strand

(all coordinates are 0-based)
'''

import os
import sys
from ngsutils.gtf import GTF
from ngsutils.support.eta import ETA


class JuncsWriter(object):
    def __init__(self):
        self.junctions = set()

    def write(self, chrom, left, right, strand):
        if not (chrom, left, right, strand) in self.junctions:
            self.junctions.add((chrom, left, right, strand))

    def dump(self):
        for chrom, left, right, strand in sorted(list(self.junctions)):
            sys.stdout.write('%s\n' % '\t'.join([str(x) for x in [chrom, left, right, strand]]))
        self.junctions = set()


def gtf_juncs(gtf_fname, only_canonical=False):
    gtf = GTF(gtf_fname)

    eta = ETA(gtf.fsize(), fileobj=gtf)
    writer = JuncsWriter()

    for gene in gtf.genes:
        for txpt in gene.transcripts:
            if len(txpt.exons) > 1000 or gene.gene_name == 'abParts':
                # skip IG hyper / Ab regions
                continue
            lastend = None
            for i, (start, end) in enumerate(txpt.exons):
                eta.print_status(extra='%s:%s %s #%s' % (gene.chrom, gene.start, gene.gene_name, i))
                if lastend:
                    if canonical:
                        writer.write(gene.chrom, lastend, start, gene.strand)
                    else:
                        for txpt2 in gene.transcripts:
                            for (s, e) in txpt2.exons:
                                if s > lastend:
                                    writer.write(gene.chrom, lastend, s, gene.strand)

                lastend = end - 1
        writer.dump()
    eta.done()


def usage(msg=None):
    if msg:
        print msg
    print __doc__
    print '''\
Usage: gtfutils juncts {opts} genes.gtf{.gz}

Arguments
  genes.txt       Gene model in GTF format

Options
  -canonical      Include only junctions that are in the GTF model
                  By default we export all possible junctions between
                  all exons within the same gene.
'''
    sys.exit(1)

if __name__ == '__main__':
    gtf = None
    canonical = False

    for arg in sys.argv[1:]:
        if arg == '-canonical':
            canonical = True
        elif gtf is None and os.path.exists(arg):
            gtf = arg

    if not gtf:
        usage()

    gtf_juncs(gtf, canonical)
