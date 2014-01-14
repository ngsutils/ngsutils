#!/usr/bin/env python
## category General
## desc Annotates genomic positions based on a GTF model
'''\
Annotates genomic positions based on a GTF model

For a given input tab-delimited text file with a valid reference and position column,
this script will add columns for the gene name and intron/exon/utr/intragenic position
within the gene (if the position maps to a gene).

Column counts start at 1.
'''

import sys
import os
from ngsutils.gtf import GTF
import ngsutils.support


def usage(msg=None):
    if msg:
        print '%s\n' % msg
    print __doc__
    print '''Usage: gtfutils annotate {options} filename.gtf{.gz} input.txt

Options:
    -ref num        Column with reference (default: 1)
    -pos num        Column with pos (1-based) (default: 2)

    -gene_id        Output gene_id (from GTF file)
    -transcript_id  Output transcript_id (from GTF file)
    -gene_name      Output gene_name (from GTF file)
    -gene_location  Output gene location (exon, intron, etc)

    -noheader       The first line is not a header (default: True)
'''
    sys.exit(1)


def gtf_annotate(gtf, infile, ref_col=1, pos_col=2, gene_name=False, gene_location=False, gene_id=False, transcript_id=False, header=True):
    numcols = 0

    for line in ngsutils.support.gzip_reader(infile):
        cols = line.strip().split('\t')
        if not numcols:
            numcols = len(cols)

        if header:
            if gene_id:
                cols.append('gene_id')
            if transcript_id:
                cols.append('transcript_id')
            if gene_name:
                cols.append('gene_name')
            if gene_location:
                cols.append('gene_location')
            header = False

        else:
            while len(cols) < numcols:
                cols.append('')

            ref = cols[ref_col]
            pos = int(cols[pos_col])

            gene_ids = []
            txpt_ids = []
            gene_names = []
            locs = []

            if ref and pos:
                for gene in gtf.find(ref, pos):
                    gene_names.append(gene.gene_name)
                    gene_ids.append(gene.gene_id)
                    for txpt in gene.transcripts:
                        txpt_ids.append(txpt.transcript_id)
                        found = False
                        for start, end in txpt.exons:
                            if start < pos < end:
                                if pos < txpt.start_codon:
                                    locs.append("5'UTR")
                                elif pos > txpt.stop_codon:
                                    locs.append("3'UTR")
                                else:
                                    locs.append('coding')
                                found = True
                                break
                        if not found:
                            locs.append('intron')

                if not locs:
                    locs = ['intergenic']

            if gene_id:
                cols.append(','.join(gene_ids) if gene_ids else '')
            if transcript_id:
                cols.append(','.join(txpt_ids) if txpt_ids else '')
            if gene_name:
                cols.append(','.join(gene_names) if gene_names else '')
            if gene_location:
                cols.append(','.join(locs) if locs else '')

        print '\t'.join(cols)


if __name__ == '__main__':
    gtffile = None
    infile = None
    ref_col = 0  # stored here as 0-based, arg set as 1-based
    pos_col = 1

    gene_id = False
    transcript_id = False
    gene_name = False
    gene_location = False

    header = True

    last = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif last == '-ref':
            ref_col = int(arg) - 1
            last = None
        elif last == '-pos':
            pos_col = int(arg) - 1
            last = None
        elif arg == '-noheader':
            header = False
        elif arg == '-gene_id':
            gene_id = True
        elif arg == '-transcript_id':
            transcript_id = True
        elif arg == '-gene_name':
            gene_name = True
        elif arg == '-gene_location':
            gene_location = True
        elif arg in ['-ref', '-pos']:
            last = arg
        elif not gtffile and os.path.exists(arg):
            gtffile = arg
        elif not infile and os.path.exists(arg):
            infile = arg
        else:
            print 'Unknown argument: %s' % arg

    if not gtffile:
        usage('Missing GTF file')
    if not infile:
        usage('Missing input file')
    if not (gene_name or gene_location or gene_id or transcript_id):
        usage('Missing outputs - nothing to annotate')

    gtf = GTF(gtffile)
    gtf_annotate(gtf, infile, ref_col, pos_col, gene_name, gene_location, gene_id, transcript_id, header)
