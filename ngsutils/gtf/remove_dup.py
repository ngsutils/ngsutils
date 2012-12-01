#!/usr/bin/env python
## category General
## desc Removes "_dup" entries from a GTF file
## experimental
'''In RefSeq annotations, a particular transcript may be assigned to multiple
locations. In this case, the transcript name is altered to be
transciptid_dup1. This program will remove all _dup? entries from a GTF file.
'''

import sys
import os
from ngsutils.support import gzip_reader


def gtf_remove_dup(fname, out=sys.stdout, quiet=False):
    if not quiet:
        sys.stderr.write('Reading GTF...\n')
    dup_count = 0
    good_count = 0
    for line in gzip_reader(fname, quiet=quiet):
        try:
            if line[0] == '#':
                out.write(line)
                continue

            chrom, source, feature, start, end, score, strand, frame, attrs = line.rstrip().split('\t')
            transcript_id = None
            for key, val in [x.split(' ') for x in [x.strip() for x in attrs.split(';')] if x]:
                if val[0] == '"' and val[-1] == '"':
                    val = val[1:-1]
                if key == 'transcript_id':
                    transcript_id = val

            if '_dup' in transcript_id:
                dup_count += 1
                continue

            good_count += 1
            out.write(line)

        except:
            import traceback
            sys.stderr.write('Error parsing line:\n%s\n' % line)
            traceback.print_exc()
            sys.exit(1)

    if not quiet:
        sys.stderr.write('Kept %s transcript/exon annotations\n' % good_count)
        sys.stderr.write('Removed %s duplicate transcript/exon annotations\n' % dup_count)

    return (good_count, dup_count)


def usage(msg=None):
    if msg:
        print msg
    print __doc__
    print 'Usage: gtfutils remove_dup filename.gtf'
    sys.exit(1)

if __name__ == '__main__':
    gtf = None

    for arg in sys.argv[1:]:
        if not gtf and (os.path.exists(arg) or arg == '-'):
            gtf = arg

    if not gtf:
        usage()

    gtf_remove_dup(gtf)
