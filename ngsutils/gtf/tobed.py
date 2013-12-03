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
    -introns  Each annotated intron
    -regions  Export constant / alternative regions (annotated spliced regions)
    -tss      Transcription start sites (unique)
    -txs      Transcription stop sites (unique)
    -tlss     Translational start sites (unique start codons)
    -tlxs     Translational stop sites (unique stop codons)
    -junc5    Splice junction 5' donor
    -junc3    Splice junction 3' acceptor

    -promoter length    Promoter region from the gene [length] upstream of TSS
                        
                        Note: Length may also be in the form "up,down", where
                        the promoter coordinates will be TSS-up -> TSS+down.

                        By default the "down" length is zero.

                        For example, for a gene that starts a chr1:1000 (+), using
                        "-promoter 200,100" would yield a BED region of:
                        chr1   800    1100

'''
    sys.exit(1)


def gtf_junc_5_tobed(gtf, out=sys.stdout):
    for gene in gtf.genes:
        sites = set()
        for i, txscr in enumerate(gene.transcripts):
            if gene.strand == '+':
                for j, (start, end) in enumerate(txscr.exons):
                    if j == len(txscr.exons) - 1:
                        continue

                    if not end in sites:
                        out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, end, end + 1, '%s.%s' % (gene.gene_name, i), 0, gene.strand]]))
                        sites.add(end)
            else:
                for j, (start, end) in enumerate(txscr.exons):
                    if j == 0:
                        continue

                    if not start in sites:
                        out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, start - 1, start, '%s.%s' % (gene.gene_name, i), 0, gene.strand]]))
                        sites.add(end)


def gtf_junc_3_tobed(gtf, out=sys.stdout):
    for gene in gtf.genes:
        sites = set()
        for i, txscr in enumerate(gene.transcripts):
            if gene.strand == '-':
                for j, (start, end) in enumerate(txscr.exons):
                    if j == len(txscr.exons) - 1:
                        continue

                    if not end in sites:
                        out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, end, end + 1, '%s.%s' % (gene.gene_name, i), 0, gene.strand]]))
                        sites.add(end)
            else:
                for j, (start, end) in enumerate(txscr.exons):
                    if j == 0:
                        continue

                    if not start in sites:
                        out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, start - 1, start, '%s.%s' % (gene.gene_name, i), 0, gene.strand]]))
                        sites.add(end)


def gtf_genes_tobed(gtf, out=sys.stdout):
    for gene in gtf.genes:
        out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, gene.start, gene.end, gene.gene_name, 0, gene.strand]]))


def gtf_promoter_tobed(gtf, promoter_up, promoter_down, out=sys.stdout):
    for gene in gtf.genes:
        sites = set()
        for i, txscr in enumerate(gene.transcripts):
            if gene.strand == '+':
                if not txscr.start in sites:
                    out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, txscr.start - promoter_up, txscr.start + promoter_down, '%s.%s' % (gene.gene_name, i), 0, gene.strand]]))
                    sites.add(txscr.start)
            else:
                if not txscr.end in sites:
                    out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, txscr.end - promoter_down, txscr.end + promoter_up, '%s.%s' % (gene.gene_name, i), 0, gene.strand]]))
                    sites.add(txscr.end)

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


def gtf_txs_tobed(gtf, out=sys.stdout):
    for gene in gtf.genes:
        sites = set()
        for i, txscr in enumerate(gene.transcripts):
            if gene.strand == '-':
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
            if not txscr.has_cds:
                continue
            if not txscr.start_codon in sites:
                out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, txscr.start_codon[0], txscr.start_codon[1], '%s.%s' % (gene.gene_name, i), 0, gene.strand]]))
                sites.add(txscr.start_codon)


def gtf_tlxs_tobed(gtf, out=sys.stdout):
    'Outputs all exons (from all transcripts)'
    for gene in gtf.genes:
        sites = set()
        for i, txscr in enumerate(gene.transcripts):
            if not txscr.has_cds:
                continue
            if not txscr.stop_codon in sites:
                out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, txscr.stop_codon[0], txscr.stop_codon[1], '%s.%s' % (gene.gene_name, i), 0, gene.strand]]))
                sites.add(txscr.start_codon)


def gtf_exons_tobed(gtf, out=sys.stdout):
    'Outputs all exons (from all transcripts)'
    for gene in gtf.genes:
        exons = set()
        for txscr in gene.transcripts:
            exons.update(txscr.exons)

        for i, (start, end) in enumerate(sorted(exons)):
            out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, start, end, '%s.e%s' % (gene.gene_name, i + 1), 0, gene.strand]]))


def gtf_introns_tobed(gtf, out=sys.stdout):
    'Outputs all introns (from all transcripts)'
    for gene in gtf.genes:
        introns = set()
        
        for txscr in gene.transcripts:
            last = None
            for start, end in txscr.exons:
                if last:
                    introns.add((last, start))
                last = end

        for i, (start, end) in enumerate(sorted(introns)):
            out.write('%s\n' % '\t'.join([str(x) for x in [gene.chrom, start, end, '%s.i.%s' % (gene.gene_name, i + 1), 0, gene.strand]]))


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
    introns = False
    regions = False
    tss = False
    tlss = False
    txs = False
    tlxs = False
    junc_5 = False
    junc_3 = False
    promoter = False
    promoter_up = 0
    promoter_down = 0
    last = None

    filename = None


    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif last == '-promoter':
            if ',' in arg:
                promoter_up, promoter_down = [int(x) for x in arg.split(',')]
            else:
                promoter_up = int(arg)
            last = None
        elif arg == '-genes':
            genes = True
        elif arg == '-exons':
            exons = True
        elif arg == '-introns':
            introns = True
        elif arg == '-regions':
            regions = True
        elif arg == '-tss':
            tss = True
        elif arg == '-tlss':
            tlss = True
        elif arg == '-txs':
            txs = True
        elif arg == '-tlxs':
            tlxs = True
        elif arg == '-junc5':
            junc_5 = True
        elif arg == '-junc3':
            junc_3 = True
        elif arg in ['-promoter']:
            promoter = True
            last = arg
        elif not filename and os.path.exists(arg):
            filename = arg

    i = 0
    for arg in [genes, exons, introns, regions, tss, tlss, txs, tlxs, junc_5, junc_3, promoter]:
        if arg:
            i += 1

    if i == 0:
        usage('You must select one [type] to export.')
    elif i > 1:
        usage('You must select *only one* [type] to export.')
    elif not filename:
        usage('Missing input file')
    elif promoter and not (promoter_down or promoter_up):
        usage('You must specify a valid promoter length!')

    gtf = GTF(filename)

    if genes:
        gtf_genes_tobed(gtf)
    elif exons:
        gtf_exons_tobed(gtf)
    elif introns:
        gtf_introns_tobed(gtf)
    elif regions:
        gtf_regions_tobed(gtf)
    elif tss:
        gtf_tss_tobed(gtf)
    elif tlss:
        gtf_tlss_tobed(gtf)
    elif txs:
        gtf_txs_tobed(gtf)
    elif tlxs:
        gtf_tlxs_tobed(gtf)
    elif junc_5:
        gtf_junc_5_tobed(gtf)
    elif junc_3:
        gtf_junc_3_tobed(gtf)
    elif promoter:
        gtf_promoter_tobed(gtf, promoter_up, promoter_down)

