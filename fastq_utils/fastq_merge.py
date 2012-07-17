#!/usr/bin/env python
## category General
## desc Merges paired FASTQ files into one file
'''
Merges two (or more) paired end FASTQ files together for combined mapping.
The files need to have the paired reads in the same order. They will be
written out as:

@name1
seq1
+
qual1
@name1
seq2
+
qual2
...

The merged file is written to stdout.
'''

import os
import sys

from fastq_utils import read_fastq


def fastq_merge(fnames, split_slashes=False):
    infiles = []

    first = True
    for fname in fnames:
        gen = read_fastq(fname, quiet=not first)
        infiles.append((fname, gen))
        first = False

    while True:
        lastname = None

        try:
            for fname, generator in infiles:
                name, seq, qual = generator.next()

                if split_slashes and '/' in name:
                    spl = name.split('/', 1)
                    name = spl[0]
                    desc = ' /%s' % spl[1]
                else:
                    cols = name.split()
                    name = cols[0]
                    if len(cols) > 1:
                        desc = cols[1]
                    else:
                        desc = ''

                if not lastname:
                    lastname = name
                elif name != lastname:
                    sys.stderr.write('Files are not paired! (error in: %s)\nExpected: %s\nGot     : %s\n' % (fname, lastname, name))
                    sys.exit(1)

                sys.stdout.write('%s %s\n%s\n+\n%s\n' % (name, desc, seq, qual))
        except:
            break


def usage():
    print __doc__
    print """Usage: fastqutils merge {-slash} file1.fastq{.gz} file2.fastq{.gz} ...

-slash    Split the read name at a '/' (Illumina paired format)
"""
    sys.exit(1)

if __name__ == '__main__':
    fnames = []
    split_slashes = False
    for arg in sys.argv[1:]:
        if arg == '-slash':
            split_slashes = True
        elif os.path.exists(arg):
            fnames.append(arg)

    if not fnames:
        usage()

    fastq_merge(fnames, split_slashes)
