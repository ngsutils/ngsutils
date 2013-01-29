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

from ngsutils.fastq import FASTQ


def generator_fetch(generators):
    while True:
        try:
            yield [gen.next() for gen in generators]
        except StopIteration:
            return


def fastq_merge(fastqs, split_slashes=False, out=sys.stdout, quiet=False):
    for reads in generator_fetch([fq.fetch(quiet=quiet if i == 0 else False) for i, fq in enumerate(fastqs)]):
        cur_name = None
        for read in reads:
            name = read.name
            comment = read.comment

            if split_slashes and '/' in name:
                spl = name.split('/', 1)
                name = spl[0]
                if read.comment:
                    comment = '/%s %s' % (spl[1], read.comment)
                else:
                    comment = '/%s' % spl[1]

            if not cur_name:
                cur_name = name
            else:
                if name != cur_name:
                    raise ValueError('Files are not paired! Expected: "%s", got "%s"!' % (cur_name, name))

            read.clone(name=name, comment=comment).write(out)


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

    if len(fnames) < 2:
        usage()

    fastqs = [FASTQ(x) for x in fnames]

    fastq_merge(fastqs, split_slashes)

    for fq in fastqs:
        fq.close()
