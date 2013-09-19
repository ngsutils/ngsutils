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
        sys.stdout.write('%s\n\n' % msg)
    sys.stdout.write(__doc__)
    sys.stdout.write('''\
Usage: bamutils pcrdup {options} infile.bam

Options:
    -frag                The reads are single-end fragments, so mark PCR
                         duplicated based only on the location of the read 
                         (not-recommended)

    -bam filename        Output BAM file with PCR duplicates marked

    -counts filename     Output of the number of reads at each position 

                         Note: this is actually the number of duplicate reads
                         at each position. If a position has multiple reads
                         mapped to it, but they are not pcr duplicates, then
                         there each will be reported separately.

    You must set either -bam or -counts (or both).

''')

    sys.exit(1)


def __flush_cur_reads(cur_reads, outbam, inbam, countfile=None):
    if cur_reads:
        for k in cur_reads:
            count = 0

            for i, (mapq, idx, r) in enumerate(sorted(cur_reads[k])[::-1]):
                count += 1
                if i > 0:
                    r.is_duplicate = True
                if outbam:
                    outbam.write(r)

            if countfile:
                countfile.write('%s\t%s\t%s\t%s\n' % (inbam.references[k[0]], k[1], k[2], count))


def pcrdup_mark(inbam, outbam, fragment=False, countfile=None):
    cur_pos = None
    cur_reads = {}

    total = 0
    unique = 0
    duplicates = 0

    dup_list = set()

    def callback(read):
        return '%s, %s, %s - %s' % (total, unique, duplicates, read.qname)

    for read in bam_iter(bamfile, callback=callback):
        if not read.is_paired or read.is_read1:
            total += 1

        if read.is_unmapped:
            __flush_cur_reads(cur_reads, outbam, inbam, countfile)
            if outbam:
                outbam.write(read)
            continue

        start_pos = (read.tid, read.pos)

        if fragment:
            dup_pos = (read.tid, read.pos, '')
        else:
            # isize is the insert length, which if this is the first read, will
            # be the right most part of the second read. If the ends of the reads
            # are trimmed for QC reasons, only the 5' pos of the first read and the 3'
            # pos of the second read will be accurate.
            
            dup_pos = (read.tid, read.pos, read.isize)

        if not cur_pos or start_pos != cur_pos:
            __flush_cur_reads(cur_reads, outbam, inbam, countfile)

            cur_pos = start_pos
            cur_reads = {}
            idx = 0

        if not fragment and (read.mate_is_unmapped or not read.is_paired or not read.is_proper_pair or read.isize < 0):
            # this is a paired file, but the mate isn't paired or proper or mapped
            # just write it out, no flags to set.

            if read.qname in dup_list:
                read.is_duplicate = True
                dup_list.remove(read.qname)

            if outbam:
                outbam.write(read)
        elif dup_pos in cur_reads:
            duplicates += 1
            if not fragment:
                dup_list.add(read.qname)
            cur_reads[dup_pos].append((read.mapq, -idx, read))
        else:
            unique += 1
            cur_reads[dup_pos] = [(read.mapq, -idx, read), ]

        idx += 1

    __flush_cur_reads(cur_reads, outbam, inbam, countfile)

    sys.stdout.write('Total reads:\t%s\n' % total)
    sys.stdout.write('Unique reads:\t%s\n' % unique)
    sys.stdout.write('PCR duplicates:\t%s\n' % duplicates)

if __name__ == '__main__':
    infile = None
    outfile = None
    countfname = None
    fragment = False

    last = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif last == '-counts':
            countfname = arg
            last = None
        elif last == '-bam':
            outfile = arg
            last = None
        elif arg in ['-counts', '-bam']:
            last = arg
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

    if not infile or not (outfile or countfname):
        usage()

    bamfile = pysam.Samfile(infile, "rb")
    bamout = None
    if outfile:
        bamout = pysam.Samfile(outfile, "wb", template=bamfile)

    if countfname:
        countfile = open(countfname, 'w')
    else:
        countfile = None

    pcrdup_mark(bamfile, bamout, fragment, countfile)

    bamfile.close()
    if bamout:
        bamout.close()
    if countfile:
        countfile.close()