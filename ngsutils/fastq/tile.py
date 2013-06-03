#!/usr/bin/env python
## category General
## desc Splits long FASTQ reads into smaller (tiled) chunks
'''
For each read in a FASTQ file, split it into smaller (overlapping) chunks.
Fragments are defined by their length and offset. For example, if the length
is 35 and the offset is 10, sub-reads will be 1->35, 11->45, 21->55, etc... If
the offset and the length are the same, then the fragments will not overlap.

Output is a set of (gzip compressed) FASTQ files
'''

import os
import sys
import gzip

from ngsutils.fastq import FASTQ


def _open_file(outbase, i, gz, quiet=False):
    if gz:
        fn = '%s.%s.fastq.gz' % (outbase, i + 1)
        tmp = os.path.join(os.path.dirname(fn), '.tmp.%s' % os.path.basename(fn))

        if not quiet:
            sys.stderr.write('Output file: %s\n' % fn)
        return (gzip.open(tmp, 'w'), tmp, fn)
    else:
        fn = '%s.%s.fastq' % (outbase, i + 1)
        tmp = os.path.join(os.path.dirname(fn), '.tmp.%s' % os.path.basename(fn))

        if not quiet:
            sys.stderr.write('Output file: %s\n' % fn)

        return (open(tmp, 'w'), tmp, fn)



def fastq_tile(fname, outbase, length, offset, gz=False, quiet=False):
    fastq = FASTQ(fname)

    outs = []
    fnames = []

    for read in fastq.fetch(quiet=quiet):
        out_idx = 0
        pos = 0
        while pos + length < len(read.seq):
            if len(outs) <= out_idx:
                fobj, tmp, fn = _open_file(outbase, out_idx, gz, quiet)
                outs.append(fobj)
                fnames.append((tmp, fn))

            read.subseq(pos, pos + length, comment="#tile:%s,%s" % (pos, pos + length)).write(outs[out_idx])
            pos += offset
            out_idx += 1

    for out in outs:
        out.close()

    fastq.close()

    for tmp, fname in fnames:
        os.rename(tmp, fname)


def usage(msg=None):
    if msg:
        print msg

    print __doc__
    print """\
Usage: fastqutils tile {opts} filename.fastq{.gz} out_template

Options:
  -len val         Length of each fragment (default: 35)
  -offset val      Offset for each fragment (default: 10)

  -gz              gzip compress the output FASTQ files

"""
    sys.exit(1)

if __name__ == '__main__':
    fname = None
    outtemplate = None
    gz = False
    length = 35
    offset = 10
    last = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        if last == '-len':
            length = int(arg)
            last = None
        elif last == '-offset':
            offset = int(arg)
            last = none
        elif arg == '-gz':
            gz = True
        elif arg in ['-len', '-offset']:
            last = arg
        elif not fname:
            if not os.path.exists(arg):
                usage("Missing file: %s" % arg)
            fname = arg
        elif not outtemplate:
            outtemplate = arg

    if not fname or not outtemplate:
        usage()

    fastq_tile(fname, outtemplate, length, offset, gz)

