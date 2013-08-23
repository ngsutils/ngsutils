#!/usr/bin/env python
## category Conversion
## desc Converts Repeatmasker.org formatted regions to FASTA
'''
Takes a repeat-masker formated output file and converts it to FASTA format.
It will optionally only output specific classes of repeats (rRNA, tRNA,
simple repeats, etc...)
'''

import sys
import os
import gzip
import pysam


def repeat2fasta(repeat_fname, ref_fname, repeat_family=None):
    if repeat_fname[-3:] == '.gz':
        repeat_f = gzip.open(repeat_fname)
    else:
        repeat_f = open(repeat_fname)
    
    ref = pysam.Fastafile(ref_fname)
    repeat_f.next()
    repeat_f.next()
    repeat_f.next()

    for line in repeat_f:
        cols = line.strip().split()
        if not cols:
            continue
        chrom = cols[4]
        start = int(cols[5]) - 1
        end = int(cols[6])
        #strand = '+' if cols[8] == '+' else '-'
        family = cols[10]
        member = cols[9]
        
        if repeat_family and family != repeat_family:
            continue

        sys.stdout.write('>%s|%s|%s:%s-%s\n' % (family, member, chrom, start, end))
        sys.stdout.write('%s\n' % wrap(ref.fetch(chrom, start, end)))

    repeat_f.close()
    ref.close()

def wrap(s, length=50):
    ar = []
    while s:
        ar.append(s[:length])
        s = s[length:]
    return '\n'.join(ar)

def usage(msg=None):
    if msg:
        sys.stdout.write('%s\n' % msg)
    sys.stdout.write('%s\n' % __doc__)
    sys.stdout.write('Usage: repeat2fasta repeatmasker.out{.gz} reference.fa [repeat-family]\n\n')
    sys.exit(1)

if __name__ == '__main__':
    repeat = None
    ref = None
    repeat_family = None

    for arg in sys.argv[1:]:
        if not repeat:
            if os.path.exists(arg):
                repeat = arg
            else:
                usage('%s missing!' % arg)
        elif not ref:
            if os.path.exists(arg):
                ref = arg
            else:
                usage('%s missing!' % arg)
        elif not repeat_family:
            repeat_family = arg

    if not repeat or not ref:
        usage()

    repeat2fasta(repeat, ref, repeat_family)
