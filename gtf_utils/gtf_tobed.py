#!/usr/bin/env python
'''Convert a GFF/GTF file to BED format

This will convert whole genes, individual exons, or expressed regions.
Expressed regions are distinct sections of exons that take into account
alternative splicing, such that each region is assigned to be 'constant' or
'alternative'.
'''

import sys
import os
from gtf_utils import GTF


def usage(msg=None):
    if msg:
        sys.stderr.write('%s\n\n' % msg)
    sys.stderr.write(__doc__)
    sys.stderr.write('''
Usage: gtfutils tobed [-genes|-exons|-regions] filename.gtf{.gz}
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


def gtf_regions_tobed(fname):
    'Outputs all regions (from all transcripts)'
    gtf = GTF(fname)
    for gene in gtf.genes:
        for i, start, end, const, names in gene.regions:
            sys.stdout.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, start, end, '%s.%s.%s' % (gene.gene_name, 'const' if const else 'alt', i), len(names), gene.strand]]))


if __name__ == '__main__':
    genes = False
    exons = False
    regions = False
    filename = None

    for arg in sys.argv[1:]:
        if arg == '-genes':
            genes = True
        elif arg == '-exons':
            exons = True
        elif arg == '-regions':
            regions = True
        elif not filename and os.path.exists(arg):
            filename = arg

    if not filename:
        usage('Missing input file')
    elif not genes and not exons and not regions:
        usage('You must select "-exons" or "-genes" or "-regions"')

    if genes:
        gtf_genes_tobed(filename)
    elif exons:
        gtf_exons_tobed(filename)
    elif regions:
        gtf_regions_tobed(filename)
