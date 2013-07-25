#!/usr/bin/env python
## category General
## desc Calculate summary statistics for a FASTQ file
'''
Calculates summary statistics for a FASTQ file.

It can calculate distribution of read lengths, average quality by base
position, if a file is in colorspace, if it contains paired end data, and
what encoding is used for the quality values (Sanger or Illumina).

Note: Any quality values less than 0 are treated as 0.
'''

import os
import sys
import collections

from ngsutils.fastq import FASTQ

StatsValues = collections.namedtuple('StatsValues', 'mean stdev min_val pct25 pct50 pct75 max_val total')


class FASTQStats(collections.namedtuple('FASTQStats', 'fastq total_reads totals lengths qualities pos_qualities')):
    @classmethod
    def _make(cls, iterable):
        result = FASTQStats(*iterable)
        result._lengthstats = None
        result._qualitystats = None
        return result

    def dump(self, out=sys.stdout, verbose=False):
        if self.fastq.is_colorspace:
            out.write("Space:\tcolorspace\n")
        else:
            out.write("Space:\tbasespace\n")

        if self.fastq.is_paired:
            out.write("Pairing:\tPaired-end (%s)\n" % self.fastq.pair_count)
        else:
            out.write("Pairing:\tFragmented\n")

        out.write("Quality scale:\t%s\n" % self.fastq.check_qualtype())
        out.write("Number of reads:\t%s\n" % self.total_reads)

        out.write('\nLength distribution\n')
        out.write('Mean:\t%s\n' % self.length_stats.mean)
        out.write('StdDev:\t%s\n' % self.length_stats.stdev)
        out.write('Min:\t%s\n' % self.length_stats.min_val)
        out.write('25 percentile:\t%s\n' % self.length_stats.pct25)
        out.write('Median:\t%s\n' % self.length_stats.pct50)
        out.write('75 percentile:\t%s\n' % self.length_stats.pct75)
        out.write('Max:\t%s\n' % self.length_stats.max_val)
        out.write('Total:\t%s\n' % self.length_stats.total)

        if verbose:
            out.write('\n')
            tmp = []
            for idx, count in enumerate(self.lengths):
                idx += 1
                if count:
                    tmp.append((idx, count))

            for idx, count in sorted(tmp)[::-1]:
                    out.write("%s\t%s\n" % (idx, count))

        out.write('\nQuality distribution\n')
        out.write('pos\tmean\tstdev\tmin\t25pct\t50pct\t75pct\tmax\tcount\n')

        for pos, stats in enumerate(self.quality_stats):
            if not stats:
                continue
            out.write('%s\t' % (pos))
            out.write('\t'.join([str(x) for x in stats]))
            out.write('\n')

        out.write('\nAverage quality string\n')

        for i, x in list(enumerate(self.qualities))[1:]:
            q = x / self.totals[i]
            if q > 0:
                out.write(chr(q + 33))
            else:
                out.write('~')

        if verbose:
            out.write('\n\nPosition\tAverage\n')
            for i, q in enumerate(self.qualities):
                out.write('%s\t%s\n' % (i, (q / self.totals[i])))

        out.write('\n')

    @property
    def length_stats(self):
        if not self._lengthstats:
            self._lengthstats = stats_counts(self.lengths)
        return self._lengthstats

    @property
    def quality_stats(self):
        if not self._qualitystats:
            self._qualitystats = []
            for pos, quals in enumerate(self.pos_qualities):
                if not quals:
                    self._qualitystats.append(None)
                else:
                    self._qualitystats.append(stats_counts(quals))
        return self._qualitystats


def fastq_stats(fastq, quiet=False):
    lengths = []    # how many reads are exactly this length?
    posquals = []   # accumulator of quality values for each position
                    # (not all the values, but an accumulator for each value at each position)
    qualities = []  # quality accumulator for each position

    total = []  # how many reads are at least this length?
                # (used for dividing an accumulator later)

    total_reads = 0
    line = 0
    try:
        for read in fastq.fetch(quiet=quiet):
            # Note: all lengths are based on the FASTQ quality score, which
            # will be the correct length for base- and color-space files. The
            # sequence may have a prefix in color-space files

            line += 4
            total_reads += 1

            while len(total) <= len(read.qual):
                total.append(0)

            for x in xrange(len(read.qual)):
                total[x + 1] += 1

            while len(lengths) <= len(read.qual):
                lengths.append(0)
                qualities.append(0)
                posquals.append([])

            lengths[len(read.qual)] += 1

            for idx, q in enumerate([ord(x) - 33 for x in read.qual]):
                if q < 0:
                    q = 0
                qualities[idx + 1] += q
                while len(posquals[idx + 1]) <= q:
                    posquals[idx + 1].append(0)
                posquals[idx + 1][q] += 1

    except KeyboardInterrupt:
        pass

    return FASTQStats._make([fastq, total_reads, total, lengths, qualities, posquals])


def stats_counts(counts):
    '''
    Takes a list of counts and calculates stats

    For example:
        [0, 1, 2, 3, 4, 5, 6] would mean:
        there were no zeros, (1) one, (2) twos, (3) threes, etc...
    '''

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
        sdacc += count * ((idx - mean) ** 2)
        acc += count
        if not pct25 and acc / total > 0.25:
            pct25 = idx
        if not pct50 and acc / total > 0.5:
            pct50 = idx
        if not pct75 and acc / total > 0.75:
            pct75 = idx

    if total > 2:
        stdev = (sdacc / (total - 1)) ** 0.5
    else:
        stdev = 0.0
    return StatsValues(mean, stdev, min_val, pct25, pct50, pct75, max_val, total)


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

    fq = FASTQ(fname)
    stats = fastq_stats(fq)
    stats.dump(verbose=verbose)
    fq.close()
