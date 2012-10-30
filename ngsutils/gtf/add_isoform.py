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
from eta import ETA
from ngsutils.support.ngs_utils import gzip_opener


def gtf_add_isoform(gtf, iso):
    isoforms = {}

    sys.stderr.write('Reading isoforms...\n')
    with gzip_opener(iso) as f:
        if iso != '-':
            eta = ETA(os.stat(iso).st_size, fileobj=f)

        for line in f:
            if line[0] == '#':
                continue
            cols = line.rstrip().split('\t')
            isoforms[cols[1]] = cols[0]
            if iso != '-':
                eta.print_status()
        if iso != '-':
            eta.done()

    sys.stderr.write('Reading/Writing GTF...\n')
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

                if transcript_id in isoforms:
                    attrs = '%s isoform_id "%s";' % (attrs, isoforms[transcript_id])

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
