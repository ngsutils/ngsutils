#!/usr/bin/env python
'''
Removes runs of homopolymers in a FASTQ file

This will remove all runs of repeated nucleotides in a FASTQ file. This can be
used to help map Ion Torrent or 454 reads, which are succeptible to errant
homopolymer calls.

This writes one new file, which is a FASTQ file with the homopolymers removed.
There is no way to recover the original sequences (unlike the FASTA based
utilities)
'''

import os
import sys
import gzip

from fastq_utils import read_fastq


def fastq_strip_homopolymer(fname, outname, gz=False):
    if outname == '-':
        out = sys.stdout
    elif gz:
        out = gzip.open(outname, 'w')
    else:
        out = open(outname, 'w')

    for name, seq, qual in read_fastq(fname):
        out.write('@%s\n' % name)
        outseq = []
        outqual = []

        if seq[0] in 'atcgATCG' and seq[1] in '01234.':
            # colorspace w/prefix
            outseq.append(seq[0])
            seq = seq[1:]

        lastbase = None
        for base, q in zip(seq, qual):
            if base != lastbase or base in '.N4':
                # if the call is a run of N's, still output all
                outseq.append(base)
                outqual.append(q)
                lastbase = base

        out.write('%s\n+%s' % (''.join(outseq), ''.join(outqual)))

    if outname != '-':
        out.close()


def usage():
    print __doc__
    print '''Usage: sequtils strip_homopolymer {-z} infile.fastq outfile.fastq

Arguments:
    infile.fastq    Input FASTQ file (can be gzipped)
    outfile.fastq   Output FASTQ file (- for stdout)

Options:
    -z              Compress the output file with gzip (not stdout)
'''
    sys.exit(1)


if __name__ == '__main__':
    fname = None
    outname = None
    gz = False
    for arg in sys.argv[1:]:
        if arg == '-z':
            gz = True
        elif not fname:
            if not os.path.exists(arg):
                print "Missing file: %s" % arg
                usage()
            fname = arg
        elif not outname:
            outname = arg

    if not fname or not outname:
        usage()

    fastq_strip_homopolymer(fname, outname, gz)
