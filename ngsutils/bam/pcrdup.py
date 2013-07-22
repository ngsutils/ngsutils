#!/usr/bin/env python
## category General
## desc Find and mark PCR duplicates
## experimental
'''
For a BAM file, find and mark all possible PCR duplicates. This is meant to be
used primarily with paired-end reads, since these have better resolution to
verify that we care seeing a legitimate PCR duplicate and not just reads that
happen to start at the same location.

The orientation for paired-end reads is assumed to be "FR" (forward-reverse).

Note: The BAM file must be sorted in order to find duplicates. For paired-end
      reads, the the proper-pair (0x4) flag must be set and the isize/tlen
      field must be correctly calculated.
'''

import sys
import os
import pysam
from ngsutils.bam import bam_iter


def usage(msg=None):
    if msg:
        sys.stderr.write('%s\n\n' % msg)
    sys.stderr.write(__doc__)
    sys.stderr.write('''\
Usage: bamutils pcrdup {options} infile.bam outfile.bam

Options:
    -frag    The reads are single-end fragments, so mark PCR duplicated
             based only on the location of the read (not-recommended)

''')

    sys.exit(1)


def __flush_cur_reads(cur_reads, outbam):
    if cur_reads:
        for k in cur_reads:
            for i, (mapq, idx, r) in enumerate(sorted(cur_reads[k])[::-1]):
                if i > 0:
                    r.is_duplicate = True
                outbam.write(r)


def pcrdup_mark(inbam, outbam, fragment=False):
    cur_pos = None
    cur_reads = {}

    total = 0
    mapped = 0
    duplicates = 0

    dup_list = set()

    for read in bam_iter(bamfile):
        total += 1

        if read.is_unmapped:
            __flush_cur_reads(cur_reads, outbam)
            outbam.write(read)
            continue

        start_pos = (read.tid, read.pos)

        if fragment:
            dup_pos = (read.tid, read.pos)
        else:
            # isize is the insert length, which if this is the first read, will
            # be the right most part of the second read. If the ends of the reads
            # are trimmed for QC reasons, only the 5' pos of the first read and the 3'
            # pos of the second read will be accurate.
            
            dup_pos = (read.tid, read.pos, read.isize)

        if not cur_pos or start_pos != cur_pos:
            __flush_cur_reads(cur_reads, outbam)

            cur_pos = start_pos
            cur_reads = {}
            idx = 0

        if not fragment and (read.mate_is_unmapped or not read.is_paired or not read.is_proper_pair or read.isize < 0):
            # this is a paired file, but the mate isn't paired or proper or mapped
            # just write it out, no flags to set.

            if read.qname in dup_list:
                read.is_duplicate = True
                dup_list.remove(read.qname)

            outbam.write(read)
        elif dup_pos in cur_reads:
            mapped += 1
            duplicates += 1
            if not fragment:
                dup_list.add(read.qname)
            cur_reads[dup_pos].append((read.mapq, -idx, read))
        else:
            mapped += 1
            cur_reads[dup_pos] = [(read.mapq, -idx, read), ]

        idx += 1

    __flush_cur_reads(cur_reads, outbam)

    sys.stdout.write('Total reads:\t%s\n' % total)
    sys.stdout.write('Proper pairs:\t%s\n' % mapped)
    sys.stdout.write('PCR duplicates:\t%s\n' % duplicates)

if __name__ == '__main__':
    infile = None
    outfile = None
    fragment = False

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif arg == '-frag':
            fragment = True
        elif not infile:
            if os.path.exists(arg):
                infile = arg
            else:
                usage("%s doesn't exist!" % arg)
        elif not outfile:
            if not os.path.exists(arg):
                outfile = arg
            else:
                usage("%s exists! Not overwriting file." % arg)

    if not infile or not outfile:
        usage()

    bamfile = pysam.Samfile(infile, "rb")
    outfile = pysam.Samfile(outfile, "wb", template=bamfile)

    pcrdup_mark(bamfile, outfile, fragment)

    bamfile.close()
    outfile.close()