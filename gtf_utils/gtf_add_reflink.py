#!/usr/bin/env python
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
from support.eta import ETA
from support.ngs_utils import  gzip_opener


def gtf_addreflink(gtf, reflink):
    link_values = {}

    sys.stderr.write('Reading refLink...\n')
    with gzip_opener(reflink) as f:
        if reflink != '-':
            eta = ETA(os.stat(reflink).st_size, fileobj=f)
        for line in f:
            cols = line.rstrip().split('\t')
            link_values[cols[2]] = (cols[0], cols[6])
            if reflink != '-':
                eta.print_status()
        if reflink != '-':
            eta.done()

    sys.stderr.write('Reading GTF...\n')
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

                if transcript_id in link_values:
                    extra = 'gene_name "%s"; isoform_id "%s";' % link_values[transcript_id]
                    attrs = '%s %s' % (attrs, extra)

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
Usage: gtfutils add_reflink filename.gtf reflink.txt
''')
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
