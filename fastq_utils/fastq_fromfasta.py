#!/usr/bin/env python
## category Conversion
## desc Converts (cs)FASTA/qual files to FASTQ format
"""
Merges a (cs)fasta file and a qual file into a FASTQ file (on stdout)

This assumes that the fasta and qual files have the same reads in the same
order (so we don't have to load all the data and merge it).
"""

import sys
import os
import gzip
from support.eta import ETA


def usage():
    print __doc__
    print """
Usage: fastqutils fromfasta {opts} filename.[cs]fasta [filename.qual]

Options:
  -q val       Use a constant value for quality (Phred char; e.g. '*' = 10)
  -tag tag     add a tag as a suffix to all read names like: readname:suffix
"""
    sys.exit(-1)


def qual_to_phred(qual):
    if type(qual) is str:
        qual = qual.strip().split()
    elif type(qual) is not type([]):
        qual = [qual, ]

    quals = []
    for q in qual:
        if not q:
            q = 0
        else:
            try:
                q = int(q)
            except ValueError, e:
                print e
                print "Line: %s" % qual
                sys.exit(2)
            if q > 93:
                q = 93
            elif q < 0:
                q = 0
        quals.append(q)

    return ''.join(['%c' % (q + 33) for q in quals])


def getline(fs):
    line = fs.readline()
    while line and line.startswith('#'):
        line = fs.readline()
    if line:
        line = line.strip()
    return line


def merge_files(fasta, qual, suffix=None, common_qual=None):
    sys.stderr.write('Merging %s and %s\n' % (os.path.basename(fasta), os.path.basename(qual) if qual else common_qual))
    if fasta.lower()[-3:] == '.gz':
        f = gzip.open(fasta)
    else:
        f = open(fasta)

    if qual:
        if qual.lower()[-3:] == '.gz':
            q = gzip.open(qual)
        else:
            q = open(qual)
    elif common_qual:
        q = None
    else:
        sys.stderr.write('Must specify a common qual value or a qual file!')
        sys.exit(-1)

    f_line = getline(f)
    if q:
        q_line = getline(q)
    else:
        q_line = None

    if ETA:
        eta = ETA(os.stat(fasta).st_size, fileobj=f)
    else:
        eta = None

    colorspace = None

    while f_line and (not q or q_line):
        if q:
            assert f_line == q_line

        name = trim_name(f_line[1:])

        if eta:
            eta.print_status(extra=name)

        if suffix:
            name = "%s:%s" % (name, suffix)
        seq = getline(f)

        if colorspace is None:
            colorspace = False
            for base in seq[1:]:
                if base in "01234":
                    colorspace = True
                    break

        if q:
            qual = qual_to_phred(getline(q))
        elif colorspace:
            qual = common_qual * (len(seq) - 1)
        else:
            qual = common_qual * len(seq)

        sys.stdout.write('@%s\n%s\n+\n%s\n' % (name, seq, qual))
        f_line = getline(f)
        if q:
            q_line = getline(q)

    if eta:
        eta.done()

    f.close()
    if q:
        q.close()


def trim_name(name):
    ''' remove trailing _F3 / _R3 (to match bfast solid2fastq script) '''
    if '_F3' in name:
        return name[:name.rindex('_F3')]
    elif '_R3' in name:
            return name[:name.rindex('_R3')]
    return name

if __name__ == '__main__':
    last = None
    tag = None
    fasta = None
    qual = None
    common_qual = None

    for arg in sys.argv[1:]:
        if last == '-tag':
            tag = arg
            last = None
        elif last == '-q':
            if len(arg) == 1:
                common_qual = arg
            else:
                sys.stderr.write("ERROR: A common qual value must be only one character\n")
                sys.exit(-1)
            last = None
        elif arg in ['-tag', '-q']:
            last = arg
        elif arg == '-h':
            usage()
        elif not fasta and os.path.exists(arg):
            fasta = arg
        elif not qual and os.path.exists(arg):
            qual = arg
        else:
            sys.stderr.write("ERROR: Unknown option: %s\n" % arg)
            usage()

    if not fasta:
        usage()

    merge_files(fasta, qual, tag, common_qual)
