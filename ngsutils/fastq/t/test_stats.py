#!/usr/bin/env python
'''
Tests for fastqutils stats
'''

import unittest
import StringIO

import ngsutils.fastq.stats
from ngsutils.fastq import FASTQ


class StatsTest(unittest.TestCase):
    def testSimple(self):
        fq = StringIO.StringIO('''\
@foo
ACGTACGT
+
;;;;;;;;
@bar
ACGTACGT
+
;;;;;;;;
@baz
ACGTACGTAC
+
;;;;;;;;;;
@quux
ACGTACGTAC
+
;;;;;;;;;;
''')
         # two reads length 8, two reads length 10
        stats = ngsutils.fastq.stats.fastq_stats(FASTQ(fileobj=fq), quiet=True)
        self.assertEqual(stats.total_reads, 4)
        self.assertEqual(stats.lengths, [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2])
        self.assertEqual(stats.totals, [0, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2])
        self.assertEqual(stats.qualities, [0, 104, 104, 104, 104, 104, 104, 104, 104, 52, 52])  # accumulator
        for q, t in zip(stats.qualities, stats.totals)[1:]:
            self.assertEqual(q / t, 26)
        self.assertEqual(stats.length_stats.mean, 9)
        self.assertEqual(stats.length_stats.min_val, 8)
        self.assertEqual(stats.length_stats.max_val, 10)

        self.assertEqual(len(stats.quality_stats), 11)

        for qvstats in stats.quality_stats:
            if not qvstats:
                continue
            self.assertEqual(qvstats.mean, 26)

    def testStatCounts(self):
        # (1) 1, (2) 2's, (3) 3's, etc...
        counts = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        stats = ngsutils.fastq.stats.stats_counts(counts)
        acc = 0
        count = 0
        for idx, val in enumerate(counts):
            idx += 1
            acc += idx * val
            count += val

        self.assertEqual(stats.min_val, 1)
        self.assertEqual(stats.max_val, 10)
        self.assertEqual(stats.pct25, 5)
        self.assertEqual(stats.pct50, 7)
        self.assertEqual(stats.pct75, 9)
        self.assertEqual(stats.total, 55)
        self.assertEqual(stats.mean, 7.0)


if __name__ == '__main__':
    unittest.main()
