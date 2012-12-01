#!/usr/bin/env python
## category General
## desc Appends isoform annotation from UCSC isoforms file
'''Adds isoform annotations to a GTF file (isoforms)

This adds isoform annotations based upon the KnownGene annotations from the
UCSC Genome Browser. Isoforms are taken from the knownIsoform table. This
table must be downloaded from the UCSC Genome Browser.

The isoforms file should be in the following format:

isoform_id {tab} transcript_id

Where all transcripts that are isoforms of the same gene will have the same
isoform_id value.

This will add the following attributes:
    isoform_id (isoform cluster id)
'''

import sys
import os
from ngsutils.support import gzip_reader


def gtf_add_isoform(gtf, iso, out=sys.stdout, quiet=False):
    isoforms = {}

    if not quiet:
        sys.stderr.write('Reading isoforms...\n')

    for line in gzip_reader(iso):
        if line[0] == '#':
            continue
        cols = line.rstrip().split('\t')
        isoforms[cols[1]] = cols[0]

    if not quiet:
        sys.stderr.write('Reading/Writing GTF...\n')

    for line in gzip_reader(gtf):
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

            if transcript_id in isoforms:
                attrs = '%s isoform_id "%s";' % (attrs, isoforms[transcript_id])

            out.write('\t'.join([chrom, source, feature, start, end, score, strand, frame, attrs]))
            if comment:
                out.write('\t%s' % comment)
            out.write('\n')
        except:
            import traceback
            sys.stderr.write('Error parsing line:\n%s\n' % line)
            traceback.print_exc()
            sys.exit(1)


def usage(msg=None):
    if msg:
        print msg
    print __doc__
    print 'Usage: gtfutils add_isoform filename.gtf knownIsoform.txt'
    sys.exit(1)

if __name__ == '__main__':
    gtf = None
    iso = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if not gtf and (os.path.exists(arg) or arg == '-'):
            gtf = arg
        elif not iso and (os.path.exists(arg) or arg == '-'):
            iso = arg

    if not gtf or not iso:
        usage()
    if gtf == '-' and iso == '-':
        usage('Both GTF and Isoform files can not be from stdin')

    gtf_add_isoform(gtf, iso)
