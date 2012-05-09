#!/usr/bin/env python
'Sorts a FASTQ file by name or sequence'


def usage(msg=None):
    if msg:
        sys.stderr.write('%s\n\n' % msg)

    sys.stderr.write(__doc__)
    sys.stderr.write("""
Usage: fastqutils sort {opts} infile.fastq{.gz} outfile.fastq{.gz}

Options:
    -name    Sort by name (default)
    -seq     Sort by sequence
    -f       Force overwriting the outfile (if it exists)
    -z       Compress the outfile
    -T dir   Temporary directory (default to current)
    -c num   The number of reads to put into each temporary file (chunk size)
             (default: 100,000)
""")
    sys.exit(1)

import sys
import os
import gzip
from fastq_utils import read_fastq
from support.eta import ETA


def fastq_sort(infile, outfile, sortby='name', compress=False, tmpdir=None, chunksize=100000):
    if not tmpdir:
        tmpdir = '.'
    tmpfiles = []
    reads = []
    count = 0
    for name, seq, qual in read_fastq(infile):
        count += 1
        if sortby == 'name':
            reads.append((name, seq, qual))
        else:
            reads.append((seq, name, qual))

        if len(reads) > chunksize:
            tmpname = os.path.join(tmpdir, '.%s.fq.%s.%s.gz' % (os.path.basename(infile), len(tmpfiles), os.getpid()))
            tmpfiles.append(tmpname)
            _write_tmp(reads, tmpname)
            reads = []

    tmpname = os.path.join(tmpdir, '.%s.fq.%s.%s.gz' % (os.path.basename(infile), len(tmpfiles), os.getpid()))
    tmpfiles.append(tmpname)
    _write_tmp(reads, tmpname)

    tmpoutname = os.path.join(os.path.dirname(outfile), '.%s.tmp' % os.path.basename(outfile))
    if compress:
        out = gzip.open(tmpoutname, 'w')
    else:
        out = open(tmpoutname, 'w')

    tmpreaders = []
    merge = []

    for t in tmpfiles:
        reader = read_fastq(t, quiet=True)
        name, seq, qual = reader.next()

        if sortby == 'name':
            merge.append((name, seq, qual, len(tmpreaders)))
        else:
            merge.append((seq, name, qual, len(tmpreaders)))

        tmpreaders.append(reader)

    eta = ETA(count)
    i = 0
    while merge:
        i += 1
        merge.sort()
        tup = merge[0]
        merge = merge[1:]
        if sortby == 'name':
            name, seq, qual, idx = tup
        else:
            seq, name, qual, idx = tup

        eta.print_status(i, extra=name)
        out.write('%s\n%s\n+\n%s\n' % (name, seq, qual))

        try:
            name, seq, qual = tmpreaders[idx].next()
            if sortby == 'name':
                merge.append((name, seq, qual, idx))
            else:
                merge.append((seq, name, qual, idx))
        except StopIteration:
            pass

    out.close()
    eta.done()
    for tmp in tmpfiles:
        os.unlink(tmp)
    os.rename(tmpoutname, outfile)


def _write_tmp(reads, tmpname):
    reads.sort()
    tmp = gzip.open(tmpname, 'w')

    for tup in reads:
        if sortby == 'name':
            name, seq, qual = tup
        else:
            seq, name, qual = tup

        tmp.write('%s\n%s\n+\n%s\n' % (name, seq, qual))
    tmp.close()

if __name__ == '__main__':
    sortby = 'name'
    compress = False
    force = False
    infile = None
    outfile = None
    tmpdir = None
    last = None
    chunksize = 100000

    for arg in sys.argv[1:]:
        if last == '-T':
            tmpdir = arg
            last = None
        elif last == '-c':
            chunksize = int(arg)
            last = None
        elif arg in ['-T', '-c']:
            last = arg
        elif arg == '-z':
            compress = True
        elif arg == '-f':
            force = True
        elif arg == '-name':
            sortby = 'name'
        elif arg == '-seq':
            sortby = 'seq'
        elif not infile and (os.path.exists(arg) or arg == '-'):
            infile = arg
        elif not outfile:
            outfile = arg
        else:
            usage('Unknown option: %s' % arg)

    if not infile:
        usage('Missing input file')
    elif not outfile:
        usage('Missing output file')
    elif not force and os.path.exists(outfile):
        usage('Will not overwrite output file without -f force')

    fastq_sort(infile, outfile, sortby, compress, tmpdir, chunksize)
