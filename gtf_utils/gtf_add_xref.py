#!/usr/bin/env python
'''Adds gene name annotations to a GTF file (xref)

This adds gene name annotations based upon the KnownGene annotations from the
UCSC Genome Browser. Gene names will be taken from the kgXref table. This
table must be downloaded separately from the UCSC Genome Browser.

This assumes that the file will be in tab-delimited format and that there is
one line for each transcript. You may specify which column represents the gene
name. In the standard "kgXref.txt" file, this is column #5.

This will add the following attributes:
    gene_name
'''

import sys
import os
from support.eta import ETA
from support.ngs_utils import  gzip_opener


def gtf_add_xref(gtf, xref, column=4):
    gene_names = {}

    sys.stderr.write('Reading xref...\n')
    with gzip_opener(xref) as f:
        if xref != '-':
            eta = ETA(os.stat(xref).st_size, fileobj=f)
        for line in f:
            if line[0] == '#':
                continue
            cols = line.rstrip().split('\t')
            gene_names[cols[0]] = cols[column]
            if xref != '-':
                eta.print_status()
        if xref != '-':
            eta.done()

    sys.stderr.write('Reading/writing GTF...\n')
    with gzip_opener(gtf) as f:
        if gtf != '-':
            eta = ETA(os.stat(gtf).st_size, fileobj=f)
        for line in f:
            try:
                comment = None
                idx = line.find('#')
                if idx > -1:
                    if idx == 0:
                        sys.stdout.write(line)
                        continue
                    comment = line[idx:]
                    line = line[:-idx]
                chrom, source, feature, start, end, score, strand, frame, attrs = line.rstrip().split('\t')
                transcript_id = None
                for key, val in [x.split(' ') for x in [x.strip() for x in attrs.split(';')] if x]:
                    if val[0] == '"' and val[-1] == '"':
                        val = val[1:-1]
                    if key == 'transcript_id':
                        transcript_id = val

                if attrs[-1] != ';':
                    attrs = '%s;' % attrs

                if transcript_id in gene_names:
                    attrs = '%s gene_name "%s";' % (attrs, gene_names[transcript_id])

                sys.stdout.write('\t'.join([chrom, source, feature, start, end, score, strand, frame, attrs]))
                if comment:
                    sys.stdout.write('\t%s' % comment)
                sys.stdout.write('\n')
                if gtf != '-':
                    eta.print_status()
            except:
                import traceback
                sys.stderr.write('Error parsing line:\n%s\n' % line)
                traceback.print_exc()
                sys.exit(1)

        if gtf != '-':
            eta.done()


def usage(msg=None):
    if msg:
        sys.stderr.write('%s\n' % msg)
    sys.stderr.write(__doc__)
    sys.stderr.write('''
Usage: gtfutils add_kg {-col num} filename.gtf kgXref.txt

Options:
  -col num    The gene name is stored in column {num} (1-based)
              (default:5)
''')
    sys.exit(1)

if __name__ == '__main__':
    gtf = None
    xref = None
    column = 4

    last = None

    for arg in sys.argv[1:]:
        if last == '-col':
            column = int(arg) - 1
            last = None
        elif not gtf and (os.path.exists(arg) or arg == '-'):
            gtf = arg
        elif not xref and (os.path.exists(arg) or arg == '-'):
            xref = arg
        elif arg in ['-col']:
            last = arg

    if not gtf or not xref:
        usage()
    if gtf == '-' and xref == '-':
        usage('Both GTF and Xref files can not be from stdin')

    gtf_add_xref(gtf, xref, column)
