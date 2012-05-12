#!/usr/bin/env python
'''Convert a GFF/GTF file to BED format (genes or exons)'''

import sys
import os
from gtf_utils import GTF


def usage(msg=None):
    if msg:
        sys.stderr.write('%s\n\n' % msg)
    sys.stderr.write(__doc__)
    sys.stderr.write('''
Usage: gtfutils tobed [-genes|-exons] filename.gtf{.gz}
''')


def gtf_genes_tobed(fname):
    gtf = GTF(fname)
    for gene in gtf.genes:
        sys.stdout.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, gene.start, gene.end, gene.gene_name, 0, gene.strand]]))


def gtf_exons_tobed(fname):
    'Outputs all exons (from all transcripts)'
    gtf = GTF(fname)
    for gene in gtf.genes:
        exons = set()
        for txscr in gene.transcripts:
            exons.update(txscr.exons)

        for i, (start, end) in enumerate(sorted(exons)):
            sys.stdout.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, start, end, '%s.%s' % (gene.gene_name, i + 1), 0, gene.strand]]))


if __name__ == '__main__':
    genes = False
    exons = False
    filename = None

    for arg in sys.argv[1:]:
        if arg == '-genes':
            genes = True
        elif arg == '-exons':
            exons = True
        elif not filename and os.path.exists(arg):
            filename = arg

    if not filename:
        usage('Missing input file')
    elif not genes and not exons:
        usage('You must select "-exons" or "-genes"')

    if genes:
        gtf_genes_tobed(filename)
    elif exons:
        gtf_exons_tobed(filename)
