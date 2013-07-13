#!/usr/bin/env python
## category General
## desc Calculate the inner mate-pair distance from two BAM files
'''
Calculate the inner mate-pair distance from two BAM files

With paired-end reads, knowing the inner mate-pair distance is key for some
downstream analyses or mapping. This value can be estimated using wet-lab
techniques, but can also be found empirically. If you map each fragment
separately to a common reference (genome or transcriptome), this command will
calculate the average distance between the fragments.

Reads where one (or both) fragment is unmapped are ignored, as are pairs that
map to different references.

(reference) ===============================================================
                 (frag 1) >>>>>>>>>                      <<<<<<<< (frag 2)
                                   |--------------------|
                                  inner mate-pair distance

This will also work if the reads are aligned in the opposite orientation:

(reference) ===============================================================
                 (frag 2) >>>>>>>>>                      <<<<<<<< (frag 1)
                                   |--------------------|
                                  inner mate-pair distance

or if they overlap:

(reference) ===============================================================
                           (frag 1) >>>>>>>>>
                                         <<<<<<<< (frag 2)
                                        |----|
                                  inner mate-pair distance (negative)

'''
import sys
import os
from ngsutils.bam import bam_iter
from ngsutils.support.stats import counts_mean_stdev
import pysam


def bam_innerdist(bam1, bam2, summaryout=None):
    iter1 = bam_iter(bam1)
    iter2 = bam_iter(bam2, quiet=True)

    distances = {}
    total = 0
    proper = 0

    orientation_count = {
        '+/-': 0,
        '-/+': 0,
        '+/+': 0,
        '-/-': 0,
    }

    read1_last = None
    read2_last = None
    read1 = None
    read2 = None

    while True:
        try:
            while not read1 or read1_last == read1.qname:
                read1 = iter1.next()
            while not read2 or read2_last == read2.qname:
                read2 = iter2.next()
        except StopIteration:
            break

        if read1.qname != read2.qname:
            raise ValueError("Error: BAM files aren't properly paired! (%s, %s)\n" % (read1.qname, read2.qname))

        read1_last = read1.qname
        read2_last = read2.qname

        total += 1

        if read1.is_unmapped or read2.is_unmapped or read1.tid != read2.tid:
            continue

        proper += 1

        if read1.pos < read2.pos:
            dist = read2.pos - read1.aend
        else:
            dist = read1.pos - read2.aend

        if summaryout:
            summaryout.write('%s\n' % dist)

        if not dist in distances:
            distances[dist] = 1
        else:
            distances[dist] += 1

        orientation = '%s/%s' % ('-' if read1.is_reverse else '+', '-' if read2.is_reverse else '+')

        orientation_count[orientation] += 1

    mean, stdev = counts_mean_stdev(distances)

    return total, proper, mean, stdev, orientation_count


def usage():  # pragma: no cover
    print __doc__
    print """\
Usage: bamutils innerdist filename1.bam filename2.bam

Options:
  -summary filename       Write all distances out to a file

Note: BAM files must be paired and they must be mapped to the
      same reference and reads must be in the same order.

"""
    sys.exit(1)

if __name__ == "__main__":  # pragma: no cover
    fname1 = None
    fname2 = None
    summary = None

    last = None
    for arg in sys.argv[1:]:
        if last == '-summary' and not os.path.exists(arg):
            summary = arg
            last = None
        elif not fname1 and os.path.exists(arg):
            fname1 = arg
        elif not fname2 and os.path.exists(arg):
            fname2 = arg
        elif arg in ['-summary']:
            last = arg
        else:
            usage()

    if not fname1 or not fname2:
        usage()

    bam1 = pysam.Samfile(fname1, "rb")
    bam2 = pysam.Samfile(fname2, "rb")

    if summary:
        fobj = open(summary, 'w')
    else:
        fobj = None

    try:
        total, proper, mean, stdev, o_count = bam_innerdist(bam1, bam2, fobj)
    except ValueError, e:
        print e
        sys.exit(1)
    finally:
        bam1.close()
        bam2.close()

        if fobj:
            fobj.close()

    print "Total reads:\t%s" % total
    print "Proper pairs:\t%s" % proper
    print ""
    print "Mean:\t%s" % mean
    print "Stdev:\t%s" % stdev

    print "Read orientation:"
    for o in o_count:
        if o_count[o]:
            print "  %s:\t%s" % (o, o_count[o])
