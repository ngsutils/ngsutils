#!/usr/bin/env python
'Index a sequence-sorted FASTQ file for prefix based searching'

import sys
import os
import support
from support.eta import ETA
from fastq_utils import is_colorspace_fastq


def usage(msg=None):
    if msg:
        sys.stderr.write('%s\n\n' % msg)

    sys.stderr.write(__doc__)
    sys.stderr.write('''
Usage: fastqutils index filename.fastq{.gz}

''')
    sys.exit(1)


def fastq_index(fname, depth=12):
    lastseq = None
    lastprefix = None
    iscolor = is_colorspace_fastq(fname)

    outname = '%s.idx' % fname
    out = open(outname, 'w')

    with open(fname) as f:
        eta = ETA(os.stat(fname).st_size, fileobj=f)

        while f:
            try:
                pos = f.tell()
                f.next()  # name
                seq = f.next().strip()
                f.next()  # +
                f.next()  # qual

                if iscolor and seq[0] in 'ATCGatgc':
                    prefix = seq[1:depth + 1]
                else:
                    prefix = seq[:depth]

                eta.print_status(extra=prefix)
                if lastseq and seq < lastseq:
                    sys.stderr.write("FASTQ file isn't sorted in sequence order! Aborting!\n")
                    out.close()
                    os.unlink(outname)
                    sys.exit(2)

                if lastprefix != prefix:
                    out.write('%s\t%s\n' % (prefix, pos))
                    lastprefix = prefix

            except KeyboardInterrupt:
                break
            except:
                pass

        eta.done()
        out.close()


if __name__ == '__main__':
    fname = None
    for arg in sys.argv[1:]:
        if not fname and os.path.exists(arg):
            fname = arg

    if not fname:
        usage('Missing input file')
    if fname[-3:] == '.gz':
        usage('Compressed FASTQ files can not be indexed')

    fastq_index(fname)
