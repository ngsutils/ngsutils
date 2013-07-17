#!/usr/bin/env python
## category Conversion
## desc Convert a GFF/GTF file to BED format
'''
Convert a GFF/GTF file to BED format

This will convert whole genes, individual exons, or expressed regions.
Expressed regions are distinct sections of exons that take into account
alternative splicing, such that each region is assigned to be 'constant' or
'alternative'.
'''

import sys
import os
from ngsutils.gtf import GTF


def usage(msg=None):
    if msg:
        print '%s\n' % msg
    print __doc__
    print '''\
Usage: gtfutils tobed [type] filename.gtf{.gz}

Where type is one of:
    -genes    The gene from start to end (including introns)
    -exons    Each annotated exon
    -regions  Export constant / alternative regions (annotated splice regions)
    -tss      Unique transcription start sites

'''
    sys.exit(1)


def gtf_genes_tobed(gtf, out=sys.stdout):
    for gene in gtf.genes:
        out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, gene.start, gene.end, gene.gene_name, 0, gene.strand]]))


def gtf_tss_tobed(gtf, out=sys.stdout):
    for gene in gtf.genes:
        if gene.strand == '+':
            out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, gene.start, gene.start + 3, gene.gene_name, 0, gene.strand]]))
        else:
            out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, gene.end - 3, gene.end, gene.gene_name, 0, gene.strand]]))


def gtf_exons_tobed(gtf, out=sys.stdout):
    'Outputs all exons (from all transcripts)'
    for gene in gtf.genes:
        exons = set()
        for txscr in gene.transcripts:
            exons.update(txscr.exons)

        for i, (start, end) in enumerate(sorted(exons)):
            out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, start, end, '%s.%s' % (gene.gene_name, i + 1), 0, gene.strand]]))


def gtf_regions_tobed(gtf, out=sys.stdout):
    'Outputs all regions (from all transcripts)'
    for gene in gtf.genes:
        for i, start, end, const, names in gene.regions:
            source_count = 0
            for n in names.split(','):
                source_count += 1
            out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, start, end, '%s.%s.%s' % (gene.gene_name, 'const' if const else 'alt', i), source_count, gene.strand]]))


if __name__ == '__main__':
    genes = False
    exons = False
    regions = False
    tss = False
    filename = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if arg == '-genes':
            genes = True
        elif arg == '-exons':
            exons = True
        elif arg == '-regions':
            regions = True
        elif arg == '-tss':
            tss = True
        elif not filename and os.path.exists(arg):
            filename = arg

    if not genes and not exons and not regions and not tss:
        usage('You must select "-exons" or "-genes" or "-regions" or "-tss"')
    elif not filename:
        usage('Missing input file')

    gtf = GTF(filename)

    if genes:
        gtf_genes_tobed(gtf)
    elif exons:
        gtf_exons_tobed(gtf)
    elif regions:
        gtf_regions_tobed(gtf)
    elif tss:
        gtf_tss_tobed(gtf)
