#!/usr/bin/env python
## category General
## desc Appends isoform/name annotation from RefSeq/refLink
'''Adds gene name and isoform annotations to a GTF file (refLink)

This adds isoform and gene name annotations based upon the refLink table from
the UCSC Genome Browser. This assumes that the GTF transcript_id is the same
as the mrnaAcc field in the refLink table. Any transcripts with the same
locuslink/gene id will be treated as isoforms of each other.

This will add the following attributes:
    gene_name
    isoform_id (NCBI LocusLinkID/GeneID)
'''

import sys
import os
from ngsutils.support import gzip_reader


def gtf_addreflink(gtf, reflink, out=sys.stdout, quiet=False):
    link_values = {}

    if not quiet:
        sys.stderr.write('Reading refLink...\n')

    for line in gzip_reader(reflink):
        cols = line.rstrip().split('\t')
        link_values[cols[2]] = (cols[0], cols[6])

    if not quiet:
        sys.stderr.write('Reading GTF...\n')

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

            if transcript_id in link_values:
                extra = 'gene_name "%s"; isoform_id "%s";' % link_values[transcript_id]
                attrs = '%s %s' % (attrs, extra)

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
    print 'Usage: gtfutils add_reflink filename.gtf reflink.txt'
    sys.exit(1)

if __name__ == '__main__':
    gtf = None
    reflink = None

    for arg in sys.argv[1:]:
        if not gtf and (os.path.exists(arg) or arg == '-'):
            gtf = arg
        elif not reflink and (os.path.exists(arg) or arg == '-'):
            reflink = arg

    if not gtf or not reflink:
        usage()
    if gtf == '-' and reflink == '-':
        usage('Both GTF and reflink files can not be from stdin')

    gtf_addreflink(gtf, reflink)
