#!/usr/bin/env python
## category General
## desc Calculate summary statistics for a FASTQ file
'''
Calculates summary statistics for a FASTQ file.

It can calculate distribution of read lengths, average quality by base
position, if a file is in colorspace, if it contains paired end data, and
what encoding is used for the quality values (Sanger or Illumina).
'''

import os
import sys

from fastq_utils import read_fastq, is_colorspace_fastq, is_paired_fastq, fastq_qualtype


def fastq_stats(fname, verbose=False):
    cs = is_colorspace_fastq(fname)
    if cs:
        print "Space:\tcolorspace"
    else:
        print "Space:\tnucleotide"

    pairs = is_paired_fastq(fname)
    if pairs > 0:
        print "Pairing:\tPaired-end (%s)" % pairs
    else:
        print "Pairing:\tFragmented"

    qual_totals = fastq_qualtype(fname)

    print "Quality scale:\t%s" % qual_totals[-1][1]
    if verbose:
        print ' '.join(['(%s,%s)' % (q[1], q[0]) for q in qual_totals])

    lengths = []
    posquals = []
    qualities = []
    total = []
    total_reads = 0
    line = 0
    try:
        for name, seq, qual in read_fastq(fname):
            if not name[0] == '@':
                print 'Invalid formatted record [line %s]' % line
                break

            if cs:
                if len(seq) != len(qual) + 1:
                    print 'Seq / qual mismatch [line %s]' % line
                    break
            else:
                if len(seq) != len(qual):
                    print 'Seq / qual mismatch [line %s]' % line
                    break

            line += 4
            total_reads += 1

            while len(total) < len(qual) + 1:
                total.append(0)
            for x in xrange(len(qual) + 1):
                total[x] += 1

            while len(qual) > len(lengths) - 1:
                lengths.append(0)
                qualities.append(0)
                posquals.append([])

            lengths[len(qual)] += 1

            for idx, q in enumerate([ord(x) for x in qual]):
                qualities[idx] += q
                while len(posquals[idx]) < (q - 32):
                    posquals[idx].append(0)
                posquals[idx][q - 33] += 1

    except KeyboardInterrupt:
        pass

    print "Number of reads:\t%s" % total_reads
    print ""

    mean, stdev, min_val, pct25, pct50, pct75, max_val = stats_counts(lengths)

    print "Length distribution"
    print 'Mean:\t%s' % mean
    print 'StdDev:\t%s' % stdev
    print 'Min:\t%s' % min_val
    print '25 percentile:\t%s' % pct25
    print 'Median:\t%s' % pct50
    print '75 percentile:\t%s' % pct75
    print 'Max:\t%s' % max_val

    if verbose:
        print ''
        for idx, count in enumerate(lengths[::-1]):
            if count:
                print "%s\t%s" % (len(lengths) - idx - 1, count)

    print "Quality distribution"
    print "pos\tmean\tstdev\tmin\t25pct\t50pct\t75pct\tmax"

    for pos, quals in enumerate(posquals):
        if not quals:
            continue
        mean, stdev, min_val, pct25, pct50, pct75, max_val = stats_counts(quals)

        sys.stdout.write('%s\t' % (pos + 1))
        sys.stdout.write('\t'.join([str(x) for x in stats_counts(quals)]))
        sys.stdout.write('\n')

    print ""
    print "Quality by position"

    for i, x in enumerate(qualities):
        q = x / total[i]
        if q > 33:
            sys.stdout.write(chr(q))
        else:
            sys.stdout.write('~')

    if verbose:
        print ''
        for i, q in enumerate(qualities):
            print '[%s] %s' % (i, (q / total[i]) - 33)

    print ''


def stats_counts(counts):
    acc = 0.0
    pct25 = None
    pct50 = None
    pct75 = None
    mean = None
    min_val = None
    max_val = None

    total = 0

    for idx, count in enumerate(counts):
        if not min_val and count:
            min_val = idx

        if count:
            max_val = idx

        acc += (idx * count)
        total += count

    mean = acc / total
    acc = 0.0
    sdacc = 0.0

    for idx, count in enumerate(counts):
        sdacc += ((idx - mean) ** 2)
        acc += count
        if not pct25 and acc / total > 0.25:
            pct25 = idx
        if not pct50 and acc / total > 0.5:
            pct50 = idx
        if not pct75 and acc / total > 0.75:
            pct75 = idx

    stdev = (sdacc / (total - 1)) ** 0.5
    return (mean, stdev, min_val, pct25, pct50, pct75, max_val)


def usage():
    print __doc__
    print "Usage: fastqutils stats filename.fastq{.gz}"
    sys.exit(1)


if __name__ == '__main__':
    fname = None
    verbose = False

    for arg in sys.argv[1:]:
        if arg == '-v':
            verbose = True
        elif arg == '-h':
            usage()
        elif os.path.exists(arg):
            fname = arg

    if not fname:
        usage()

    fastq_stats(fname, verbose)
