#!/usr/bin/env python
## category General
## desc Given two separately mapped paired files, re-pair the files
"""
Given two separately mapped paired-end files, re-pair the files, selecting
the most likely pairing partners based upon strand, insert distance, and
maximizing alignment scores.

It is very important that the files are either in the same order with each
read present in both files or sorted in name order.

The value of the attribute/tag given will be used to determine which reads
should be kept and which should be discarded. The tag should be a numeric
(int/float) type. More than one tag can be used. The default is 'AS+, NM-'.

The BEST pair will be kept that maximizes the tag values and otherwise
satisfies strand and distance values.
"""

import os
import sys
import pysam
import ngsutils.bam


def usage(msg=None):
    if msg:
        print msg
    print __doc__
    print """
Usage: bamutils pair {opts} out.bam read1.bam read2.bam 

Options
  -tag VAL            Tag to use to determine from which file reads will be
                      taken. (must be type :i or :f) You may have more than
                      one of these, in which case they will be sorted in
                      order. You can add a +/- at the end of the name to 
                      signify sort order (asc/desc). 
                      Default: AS+, NM-

  -fail1 fname.bam    Write all failed mappings from read1 to this file
  -fail2 fname.bam    Write all failed mappings from read1 to this file
                      (Note: -fail1 and -fail2 can be the same file.)

  -size low-high      The minimum/maximum insert size to accept. By default,
                      this will attempt to minimize the distance between
                      reads, upto the lower-bound. Any pair over the upper
                      bound will be discarded. Note: for RNA, you should
                      select a very large upper bound to account for splicing
                      (for example: 50-1000000).
                      Default: 50-1000.

"""
    sys.exit(1)


def is_valid_pair(read1, read2):
    if read1.is_unmapped or read2.is_unmapped:
        # both must be mapped
        return False

    if read1.tid != read2.tid:
        # to the same chromosome/reference
        return False

    if read1.is_reverse == read2.is_reverse:
        # in opposite orientations
        return False

    # sequenced towards each other
    if read1.pos < read2.pos and read1.is_reverse:
        return False

    if read2.pos < read1.pos and read2.is_reverse:
        return False

    return True


def find_pairs(reads1, reads2, min_size, max_size, tags):
    '''
    returns pairs, fail1, fail2
    '''

    possible = []
    fail1 = []
    fail2 = []

    valid = set()

    for r1 in reads1:
        for r2 in reads2:
            if is_valid_pair(r1, r2):
                if r1.pos < r2.pos:
                    size = r2.aend - r1.pos
                else:
                    size = r1.aend - r2.pos

                if size < min_size or size > max_size:
                    continue

                tag_val = []
                for tag in tags:
                    val = float(r1.opt(tag[:2]))

                    val += float(r2.opt(tag[:2]))
                    if tag[-1] == '+':
                        # we will sort ascending to minimize size, so + tags (AS) need to be reversed
                        val = -val
                    tag_val.append(val)

                possible.append((tag_val, size, r1, r2))
                valid.add((1, r1.tid, r1.pos))
                valid.add((2, r2.tid, r2.pos))

    for r1 in reads1:
        if not (1, r1.tid, r1.pos) in valid:
            fail1.append(r1)

    for r2 in reads2:
        if not (2, r2.tid, r2.pos) in valid:
            fail2.append(r2)

    return possible, fail1, fail2


def bam_pair(out_fname, read1_fname, read2_fname, tags=['AS+', 'NM-'], min_size=50, max_size=1000, fail1_fname=None, fail2_fname=None, quiet=False):
    bam1 = pysam.Samfile(read1_fname, "rb")
    bam2 = pysam.Samfile(read2_fname, "rb")
    out = pysam.Samfile('%s.tmp' % out_fname, "wb", template=bam1)

    fail1 = None
    fail2 = None

    if fail1_fname:
        fail1 = pysam.Samfile('%s.tmp' % fail1_fname, "wb", template=bam1)

    if fail2_fname:
        if fail2_fname == fail1_fname:
            fail2 = fail1
        else:
            fail2 = pysam.Samfile('%s.tmp' % fail2_fname, "wb", template=bam1)

    gen1 = ngsutils.bam.bam_batch_reads(bam1, quiet=quiet)
    gen2 = ngsutils.bam.bam_batch_reads(bam2, quiet=True)

    reads1 = None
    reads2 = None

    while True:
        try:
            if not reads1:
                reads1 = gen1.next()

            if not reads2:
                reads2 = gen2.next()
        except StopIteration:
            break

        if reads1[0].qname != reads2[0].qname:
            if reads1[0].qname < reads2[0].qname:
                reads1 = None
            else:
                reads2 = None
            continue

        pairs, failed_reads1, failed_reads2 = find_pairs(reads1, reads2, min_size, max_size, tags)
        if pairs:
            pairs.sort()  # default: max AS, min NM, min size

            tag_val, size, r1, r2 = pairs[0]
            best_val = (tag_val, size)
            best_pairs = []

            for tag_val, size, r1, r2 in pairs:
                if (tag_val, size) == best_val:
                    best_pairs.append((size, r1, r2))

            for size, r1, r2 in best_pairs:
                # good match! set the flags and write them out
                r1.is_paired = True
                r2.is_paired = True

                r1.is_proper_pair = True
                r2.is_proper_pair = True

                r1.is_read1 = True
                r2.is_read2 = True

                if r1.pos < r2.pos:
                    r1.tlen = size
                    r2.tlen = -size
                else:
                    r1.tlen = -size
                    r2.tlen = size

                r1.mate_is_reverse = r2.is_reverse
                r2.mate_is_reverse = r1.is_reverse

                r1.mate_is_unmapped = False
                r2.mate_is_unmapped = False

                r1.rnext = r2.tid
                r2.rnext = r1.tid

                r1.pnext = r2.pos
                r2.pnext = r1.pos

                r1tags = []
                for k,v in r1.tags:
                    if k == 'NH':
                        v = len(best_pairs)
                    r1tags.append((k,v))

                r2tags = []
                for k,v in r2.tags:
                    if k == 'NM':
                        v = len(best_pairs)
                    r2tags.append((k,v))

                r1.tags = r1tags
                r2.tags = r2tags

                out.write(r1)
                out.write(r2)

        if failed_reads1:
            for r in failed_reads1:
                r.is_paired = True
                r.is_proper_pair = False
                r.is_read1 = True
                fail1.write(r)
        if failed_reads2:
            for r in failed_reads2:
                r.is_paired = True
                r.is_proper_pair = False
                r.is_read2 = True
                fail2.write(r)

        reads1 = None
        reads2 = None


    bam1.close()
    bam2.close()

    out.close()
    os.rename('%s.tmp' % out_fname, out_fname)

    if fail1:
        fail1.close()
        os.rename('%s.tmp' % fail1_fname, fail1_fname)

    if fail2:
        if fail2_fname != fail1_fname:
            fail2.close()
            os.rename('%s.tmp' % fail2_fname, fail2_fname)

if __name__ == '__main__':
    out_fname = None
    read1_fname = None
    read2_fname = None
    fail1_fname = None
    fail2_fname = None
    min_size = 50
    max_size = 1000
    tags = []

    last = None

    for arg in sys.argv[1:]:
        if arg == '-h':
            usage()
        elif last == '-fail1':
            fail1_fname = arg
            last = None
        elif last == '-fail2':
            fail2_fname = arg
            last = None
        elif last == '-size':
            min_size, max_size = [int(x) for x in arg.split('-')]
            last = None
        elif last == '-tag':
            tags.append(arg)
            last = None
        elif arg in ['-tag', '-fail1', '-fail2', '-size']:
            last = arg
        elif not out_fname:
            out_fname = arg
        elif not read1_fname and os.path.exists(arg):
            read1_fname = arg
        elif not read2_fname and os.path.exists(arg):
            read2_fname = arg
        else:
            usage('Unknown option: %s' % arg)

    if not tags:
        tags = ['AS+', 'NM-']

    if not read1_fname or not read2_fname or not out_fname:
        usage()
    else:
        bam_pair(out_fname, read1_fname, read2_fname, tags, min_size, max_size, fail1_fname, fail2_fname)
