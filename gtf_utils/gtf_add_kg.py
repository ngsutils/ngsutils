#!/usr/bin/env python
'''Adds KnownGene annotations to a GTF file

This adds isoform and gene name annotations based upon the KnownGene
annotations from the UCSC Genome Browser. Gene names will be taken from the
kgXref table and isoforms from the knownIsoform table. These tables must be
downloaded separately from the UCSC Genome Browser.

This will add the following attributes:
    gene_name
    isoform_id (isoform cluster id)
'''

import sys
import os
from support.eta import ETA
from support.ngs_utils import  gzip_opener


def gtf_add_kg(gtf, xref, iso):
    gene_names = {}
    isoforms = {}

    sys.stderr.write('Reading xref...\n')
    with gzip_opener(xref) as f:
        eta = ETA(os.stat(xref).st_size, fileobj=f)
        for line in f:
            cols = line.rstrip().split('\t')
            gene_names[cols[0]] = cols[4]
            eta.print_status()
        eta.done()

    sys.stderr.write('Reading isoforms...\n')
    with gzip_opener(iso) as f:
        eta = ETA(os.stat(iso).st_size, fileobj=f)
        for line in f:
            cols = line.rstrip().split('\t')
            isoforms[cols[1]] = cols[0]
            eta.print_status()
        eta.done()

    sys.stderr.write('Reading GTF...\n')
    with gzip_opener(gtf) as f:
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

                extra = ''

                if transcript_id in gene_names:
                    extra = 'gene_name "%s";' % gene_names[transcript_id]

                if transcript_id in isoforms:
                    if extra:
                        extra = '%s ' % extra

                    extra = '%sisoform_id "%s";' % (extra, isoforms[transcript_id])

                if extra:
                    attrs = '%s %s' % (attrs, extra)

                sys.stdout.write('\t'.join([chrom, source, feature, start, end, score, strand, frame, attrs]))
                if comment:
                    sys.stdout.write('\t%s' % comment)
                sys.stdout.write('\n')
                eta.print_status()
            except:
                import traceback
                sys.stderr.write('Error parsing line:\n%s\n' % line)
                traceback.print_exc()
                sys.exit(1)

        eta.done()


def usage(msg=None):
    if msg:
        sys.stderr.write('%s\n' % msg)
    sys.stderr.write(__doc__)
    sys.stderr.write('''
Usage: gtfutils add_kg filename.gtf kgXref.txt knownIsoform.txt
''')
    sys.exit(1)

if __name__ == '__main__':
    gtf = None
    xref = None
    iso = None

    for arg in sys.argv[1:]:
        if not gtf and os.path.exists(arg):
            gtf = arg
        elif not xref and os.path.exists(arg):
            xref = arg
        elif not iso and os.path.exists(arg):
            iso = arg

    if not gtf or not xref or not iso:
        usage()

    gtf_add_kg(gtf, xref, iso)
