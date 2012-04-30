#!/usr/bin/env python
'''
Removes runs of homopolymers in a FASTA file

This will remove all runs of repeated nucleotides in a FASTA file. This can be
used to help map Ion Torrent or 454 reads, which are succeptible to errant
homopolymer calls.

This produces two files:
    1) a new FASTA file with the homopolymers removed
    2) an index file that allows for converting the new coordinates to the
       original genomic positions

The index is currently a gzipped text file in the format: (subject to change)
>reference
pos\trepeat_count\ttotal_offset
pos\trepeat_count\ttotal_offset
'''

import os
import sys
import gzip

from support.eta import ETA
from support.homopolymers import FASTAWriter, HPSIndex


def read_fasta_bases(fname):
    if fname[-3:] == '.gz':
        fobj = gzip.open(fname)
    else:
        fobj = open(fname)

    eta = ETA(os.stat(fname).st_size, fileobj=fobj)

    ref = None
    for line in fobj:
        eta.print_status(extra=ref)
        line = line.strip()
        if line[0] == '>':
            ref = line[1:]
        else:
            for base in line:
                yield (ref, base.upper())

    fobj.close()
    eta.done()


def seq_strip_homopolymer(fname, outfa_name=None, outidx_name=None, outtxt_name=None, suffix=None):
    if outfa_name:
        outwriter = FASTAWriter(open(outfa_name, 'w', buffering=4 * 1024 * 1024))

    if outidx_name:
        outidx = HPSIndex(outidx_name, 'w')

    if outtxt_name:
        outtxt = open(outtxt_name, 'w', buffering=4 * 1024 * 1024)

    lastref = None
    lastbase = None

    strip_pos = 0
    repeat_count = 0

    for ref, base in read_fasta_bases(fname):
        if ref != lastref:
            if outfa_name:
                if suffix:
                    outwriter.write_ref('%s%s' % (ref, suffix))
                else:
                    outwriter.write_ref(ref)

            if outidx_name:
                if suffix:
                    outidx.write_ref('%s%s' % (ref, suffix))
                else:
                    outidx.write_ref(ref)

            lastref = ref
            lastbase = None

            strip_pos = 0
            repeat_count = 0

        if base == lastbase:
            repeat_count += 1
        else:
            lastbase = base
            if repeat_count > 0:
                if outidx_name:
                    outidx.write(strip_pos, repeat_count)
                if outtxt_name:
                    if suffix:
                        outtxt.write('%s%s\t%s\t%s\n' % (ref, suffix, strip_pos, repeat_count))
                    else:
                        outtxt.write('%s\t%s\t%s\n' % (ref, strip_pos, repeat_count))

            strip_pos += 1
            repeat_count = 0
            if outfa_name:
                outwriter.write(base)

    # if this ends in a repeat...
    if repeat_count > 0:
        if outidx_name:
            outidx.write(strip_pos, repeat_count)
        if outtxt_name:
            outtxt.write('%s\t%s\t%s\n' % (ref, strip_pos, repeat_count))

    if outfa_name:
        outwriter.close()
    if outidx_name:
        outidx.close()
    if outtxt_name:
        outtxt.close()


def usage():
    print __doc__
    print """Usage: sequtils strip_homopolymer {opts} infile.fa

Options:
    ** You must select at least one of these **
    -fa  fname   Output a FASTA file with the homopolymers removed
    -idx fname   Output an index file (binary) with the location of the homopolymers
    -txt fname   Output an index file (text) with the location of the homopolymers

    -suf val     Suffix for reference names (include space for comment)
"""
    sys.exit(1)


if __name__ == '__main__':
    fname = None
    outfa = None
    outidx = None
    outtxt = None
    suffix = None
    last = None
    for arg in sys.argv[1:]:
        if last == '-suf':
            suffix = arg
            last = None
        elif last == '-fa':
            outfa = arg
            last = None
        elif last == '-idx':
            outidx = arg
            last = None
        elif last == '-txt':
            outtxt = arg
            last = None
        elif arg in ['-suf', '-fa', '-idx', '-txt']:
            last = arg
        elif not fname:
            if not os.path.exists(arg):
                print "Missing file: %s" % arg
                usage()
            fname = arg

    if not fname or (not outfa and not outidx and not outtxt):
        usage()

    seq_strip_homopolymer(fname, outfa, outidx, outtxt, suffix)
