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
    -tss      Transcription start sites (unique)
    -tlss     Translational start sites (unique start codons)

'''
    sys.exit(1)


def gtf_genes_tobed(gtf, out=sys.stdout):
    for gene in gtf.genes:
        out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, gene.start, gene.end, gene.gene_name, 0, gene.strand]]))


def gtf_tss_tobed(gtf, out=sys.stdout):
    for gene in gtf.genes:
        sites = set()
        for i, txscr in enumerate(gene.transcripts):
            if gene.strand == '+':
                if not txscr.start in sites:
                    out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, txscr.start, txscr.start + 3, '%s.%s' % (gene.gene_name, i), 0, gene.strand]]))
                    sites.add(txscr.start)
            else:
                if not txscr.end in sites:
                    out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, txscr.end - 3, txscr.end, '%s.%s' % (gene.gene_name, i), 0, gene.strand]]))
                    sites.add(txscr.end)


def gtf_tlss_tobed(gtf, out=sys.stdout):
    'Outputs all exons (from all transcripts)'
    for gene in gtf.genes:
        sites = set()
        for i, txscr in enumerate(gene.transcripts):
            if not txscr.start_codon in sites:
                out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, txscr.start_codon[0], txscr.start_codon[1], '%s.%s' % (gene.gene_name, i), 0, gene.strand]]))
                sites.add(txscr.start_codon)


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
    tlss = False
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
        elif arg == '-tlss':
            tlss = True
        elif not filename and os.path.exists(arg):
            filename = arg

    i = 0
    for arg in [genes, exons, regions, tss, tlss]:
        if arg:
            i += 1

    if i == 0:
        usage('You must select one [type] to export.')
    elif i > 1:
        usage('You must select *only one* [type] to export.')
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
    elif tlss:
        gtf_tlss_tobed(gtf)
