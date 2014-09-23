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
import gzip

from ngsutils.fastq import FASTQ, fastq_read_file
from eta import ETA


def _write_tmp(chunk, tmpdir, tmpprefix='.tmp', nogz=False):
    chunk.sort()
    tmp = tempfile.NamedTemporaryFile(prefix=tmpprefix, dir=tmpdir, delete=False)
    tmp_fname = tmp.name
    
    if not nogz: 
        fobj = gzip.GzipFile(fileobj=tmp)

        for sorter, read in chunk:
            read.write(fobj)

        fobj.close()
        tmp.close()
    else:
        for sorter, read in chunk:
            read.write(tmp)
        tmp.close()

    return tmp_fname


def fastq_sort(fastq, bysequence=False, tmpdir=None, tmpprefix='.tmp', chunksize=100000, nogz=False, out=sys.stdout, quiet=False):
    tmpfiles = []
    chunk = []

    sys.stderr.write('Sorting FASTQ file into chunks...\n')
    count = 0
    for read in fastq.fetch(quiet):
        count += 1 
        if bysequence:
            chunk.append((read.seq, read))
        else:
            chunk.append((read.name, read))

        if len(chunk) >= chunksize:
            tmpfiles.append(_write_tmp(chunk, tmpdir, tmpprefix, nogz))
            chunk = []

    if chunk:
        tmpfiles.append(_write_tmp(chunk, tmpdir, tmpprefix, nogz))

    sys.stderr.write('\nMerging chunks...\n')
    sys.stderr.flush()
    buf = [None, ] * len(tmpfiles)
    skip = [False, ] * len(tmpfiles)

    eta = ETA(count)

    j=0
    writing = True

    if nogz:
        tmpfobjs = [open(x) for x in tmpfiles]
    else:
        tmpfobjs = [gzip.open(x) for x in tmpfiles]

    while writing:
        j+=1
        eta.print_status(j)
        for i, fobj in enumerate(tmpfobjs):
            if not buf[i] and not skip[i]:
                try:
                    read = fastq_read_file(fobj)
                    if bysequence:
                        buf[i] = (read.seq, i, read)
                    else:
                        buf[i] = (read.name, i, read)
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

    for fobj in tmpfobjs:
        fobj.close()

    for tmpfile in tmpfiles:
        os.unlink(tmpfile)

def usage():
    print __doc__
    print '''
Usage: fastqutils sort {opts} filename.fastq

Options:
    -seq      Sort by read sequence (by default it sorts by name)

    -T dir    Use this directory for temporary output
    -cs num   Output this many reads in each temporary file (default: 1000000)
    -nogz     Don't compress temporary files
'''
    sys.exit(1)

if __name__ == '__main__':
    bysequence = False
    nogz = False
    tmpdir = None
    chunksize = 100000
    fname = None
    last = None
    for arg in sys.argv[1:]:
        if last == '-T':
            tmpdir = arg
            last = None
        elif last == '-cs':
            chunksize = int(arg)
        elif arg == '-nogz':
            nogz = True
        elif arg == '-seq':
            bysequence = True
        elif arg in ['-T', '-cs']:
            last = arg
        elif not fname and (os.path.exists(arg) or arg == '-'):
            fname = arg

    if not fname:
        usage()

    if not tmpdir:
        tmpdir = os.path.dirname(fname)

    fq = FASTQ(fname)
    fastq_sort(fq, bysequence=bysequence, tmpdir=tmpdir, tmpprefix='.tmp.%s' % os.path.basename(fname), chunksize=chunksize, nogz=nogz)
    fq.close()
