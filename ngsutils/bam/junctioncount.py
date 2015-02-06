#!/usr/bin/env python
## category General
## desc Counts the number of reads spanning individual junctions.
'''
Counts the number of reads that span each junction found in the BAM file.

You can specify a particular genome range to scan (like a gene region).
'''

import sys
import os
from ngsutils.bam import bam_iter, bam_open

def bam_junction_count(bam, ref=None, start=None, end=None, out=sys.stdout, quiet=False):
    last_tid = None
    junctions = {}
    for read in bam_iter(bam, ref=ref, start=start, end=end, quiet=quiet):
        if read.is_unmapped:
            continue

        if read.tid != last_tid and junctions:
            for junction in junctions:
                sys.stdout.write('%s\t%s\n' % (junction, len(junctions[junction])))
            junctions = {}
            last_tid = read.tid

        hasgap = False
        pos = read.pos
        end = None
        for op, size in read.cigar:
            if op == 0:
                pos += size
            elif op == 1:
                pass
            elif op == 2:
                pos += size
            elif op == 3:
                hasgap = True
                end = pos + size
                break
            elif op == 4:
                pos += size


        if not hasgap:
            continue

        junction = '%s:%s-%s' % (bam.references[read.tid], pos, end)
        if not junction in junctions:
            junctions[junction] = set()

        junctions[junction].add(read.qname)

    for junction in junctions:
        sys.stdout.write('%s\t%s\n' % (junction, len(junctions[junction])))


def usage(msg=""):
    if msg:
        print msg
        print
    print __doc__
    print """\
Usage: bamutils junctioncount {opts} bamfile {region}

Region should be: chr:start-end (start 1-based)

"""
    sys.exit(1)

if __name__ == "__main__":
    fname = None
    ref = None
    start = None
    end = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif not fname:
            if os.path.exists(arg):
                fname = arg
            else:
                usage("%s doesn't exist!")
        else:
            chrom, se = arg.split(':')
            start, end = [int(x) for x in se.split('-')]
            start = start - 1

    if not fname:
        usage()

    bamfile = bam_open(fname)
    bam_junction_count(bamfile, ref, start, end)
    bamfile.close()
