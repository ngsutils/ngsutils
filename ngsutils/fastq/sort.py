#!/usr/bin/env python
## category General
## desc Sorts a FASTQ file by name or sequence
'''
Sort a FASTQ file by name or sequence.

This sorts a FASTQ file into a number of smaller chunks. These chunks are then merged
together into one output written to stdout. Chunks are written to the same directory
as the original file (unless otherwise specified).
'''

import os
import sys
import tempfile

from ngsutils.fastq import FASTQ, fastq_read_file
from eta import ETA


def _write_tmp(chunk):
    chunk.sort()
    tmp = tempfile.TemporaryFile(dir=tmpdir)
    for sorter, read in chunk:
        read.write(tmp)
    tmp.seek(0)
    return tmp


def fastq_sort(fastq, byname=True, bysequence=False, tmpdir=None, chunksize=100000, out=sys.stdout, quiet=False):
    tmpfiles = []

    chunk = []
    sys.stderr.write('Sorting FASTQ file into chunks...\n')
    count = 0
    for read in fastq.fetch(quiet):
        count += 1 
        if byname:
            chunk.append((read.name, read))
        if bysequence:
            chunk.append((read.seq, read))

        if len(chunk) >= chunksize:
            tmpfiles.append(_write_tmp(chunk))
            chunk = []

    if chunk:
        tmpfiles.append(_write_tmp(chunk))

    sys.stderr.write('Merging chunks...\n')
    buf = [None, ] * len(tmpfiles)
    skip = [False, ] * len(tmpfiles)

    eta = ETA(count)

    j=0
    writing = True

    while writing:
        j+=1
        eta.print_status(j)
        for i, fobj in enumerate(tmpfiles):
            if not buf[i] and not skip[i]:
                try:
                    read = fastq_read_file(fobj)
                    if byname:
                        buf[i] = (read.name, i, read)
                    if bysequence:
                        buf[i] = (read.seq, i, read)
                except:
                    buf[i] = None
                    skip[i] = True
        
        sorted_list = buf[:]
        sorted_list.sort()
        writing = False

        for tup in sorted_list:
            if not tup:
                continue

            sorter, i, read = tup
            read.write(out)
            buf[i] = None
            writing = True
            break
    eta.done()





def usage():
    print __doc__
    print "fastqutils sort [-name | -seq] {-T dir} filename.fastq"
    sys.exit(1)

if __name__ == '__main__':
    byname = False
    bysequence = False
    tmpdir = None
    fname = None
    last = None
    for arg in sys.argv[1:]:
        if last == '-T':
            tmpdir = arg
            last = None
        elif arg in ['-T']:
            last = arg
        elif arg == '-name':
            byname = True
        elif arg == '-seq':
            bysequence = True
        elif not fname and (os.path.exists(arg) or arg == '-'):
            fname = arg

    if not fname or (not byname and not bysequence) or (byname and bysequence):
        usage()

    if not tmpdir:
        tmpdir = os.path.dirname(fname)

    fq = FASTQ(fname)
    fastq_sort(fq, byname=byname, bysequence=bysequence, tmpdir=tmpdir)
    fq.close()
