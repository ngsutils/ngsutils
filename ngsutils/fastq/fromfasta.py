#!/usr/bin/env python
## category Conversion
## desc Converts (cs)FASTA/qual files to FASTQ format
"""
Merges a (cs)fasta file and a qual file into a FASTQ file (on stdout)

This assumes that the fasta and qual files have the same reads in the same
order (so we don't have to load all the data and merge it). If there is a
constant quality value, this will automatically determine if a FASTA file is
in base-space or color-space. If the FASTA file is in color-space, the read
sequence will be examined to determine if there is a base-space prefix and the
length of the quality string adjusted accordingly.

Any '_F3' or '_R3' suffixes at the end of the read names will be removed.
"""

import sys
import os
import itertools
from ngsutils.support import FASTA


def usage():
    print __doc__
    print """
Usage: fastqutils fromfasta {opts} filename.[cs]fasta [filename.qual]

Options:
  -q val       Use a constant value for quality (Phred char; e.g. '*' = 10)
  -tag tag     add a tag as a suffix to all read names

                 Example:
                    -tag '_foo'

                 Yields:
                    @read1_foo
                    @read2_foo
                    etc...
"""
    sys.exit(-1)


def qual_to_phred(qual):
    vals = []
    for q in qual:
        if not q:
            q = 0
        else:
            try:
                q = int(q)
            except:
                raise ValueError('Error converting %s to a quality value (Line: %s)' % (q, qual))
            if q > 93:
                q = 93
            elif q < 0:
                q = 0
        vals.append(q)

    return ''.join(['%c' % (q + 33) for q in vals])


def merge_files(fasta, qual, suffix=None, common_qual=None, out=sys.stdout, quiet=False):
    colorspace = None

    if qual is None:
        qual = NullFile()

    for frec, qrec in itertools.izip(fasta.fetch(quiet=quiet), qual.fetch(quiet=False)):
        if qrec:
            if frec.name != qrec.name:
                raise ValueError("Mismatched names in FASTA and Qual files! (%s, %s)" % (frec.name, qrec.name))

        name = trim_name(frec.name)

        if suffix:
            name = "%s%s" % (name, suffix)

        if common_qual and colorspace is None:
            colorspace = False
            for base in frec.seq[1:]:
                if base in "01234":
                    colorspace = True
                    break

        if qrec:
            qual = qual_to_phred(qrec.seq.strip().split())
        elif colorspace:
            qual = common_qual * (len(frec.seq) - 1)
        else:
            qual = common_qual * len(frec.seq)

        out.write('@%s\n%s\n+\n%s\n' % (name, frec.seq, qual))


def trim_name(name):
    '''
    remove trailing _F3 / _R3 (to match bfast solid2fastq script)

    >>> trim_name('foo_F3')
    'foo'

    >>> trim_name('foo_R3')
    'foo'

    '''
    if '_F3' in name:
        return name[:name.rindex('_F3')]
    elif '_R3' in name:
            return name[:name.rindex('_R3')]
    return name


class NullFile(object):
    '''
    This object will act like a FASTA, but keep returning 'None'. This
    lets this object be used with 'zip', just like a normal FASTA iterator
    '''

    def fetch(self, *args, **kwargs):
        yield None

    def close(self):
        pass

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

    if not qual and not common_qual:
        sys.stderr.write('You must specify a quality input file or common qual value')
        usage()

    sys.stderr.write('Merging %s and %s\n' % (os.path.basename(fasta), os.path.basename(qual) if qual else common_qual))

    f = FASTA(fasta)
    q = FASTA(qual, qual=True) if qual else None

    merge_files(f, q, tag, common_qual)
    f.close()
    if q:
        q.close()
