#!/usr/bin/env python
'''
Tests for bamutils stats
'''

import os
import unittest

import ngsutils.bam
import ngsutils.bam.stats


class StatsTest(unittest.TestCase):
    def setUp(self):
        self.bam = ngsutils.bam.bam_open(os.path.join(os.path.dirname(__file__), 'test.bam'))

    def tearDown(self):
        self.bam.close()

    def testStats(self):
        stats = ngsutils.bam.stats.BamStats(self.bam)
        self.assertEqual(7, stats.total)
        self.assertEqual(6, stats.mapped)
        self.assertEqual(1, stats.unmapped)
        self.assertEqual(1, stats.flag_counts.counts[0x10])  # one reverse
        self.assertEqual(1, stats.flag_counts.counts[0x4])   # one unmapped
        self.assertEqual(0, stats.flag_counts.counts[0x2])   # zero all aligned (not paired)
        self.assertTrue('chr1' in stats.refs)   # 6 on chr1
        self.assertTrue('chr3' not in stats.refs)

    def testStatsGTF(self):
        # Add a test with a mock GTF file
        pass

if __name__ == '__main__':
    unittest.main()
