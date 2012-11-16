#!/usr/bin/env python
'''
Tests for bedutils stats
'''

import unittest
import StringIO

import ngsutils.bed.stats
from ngsutils.bed import BedFile

bedtest = BedFile(fileobj=StringIO.StringIO('''\
chr1|100|150|foo1|10|+
chr1|140|200|foo2|10|+
chr1|500|550|foo3|10|+
chr1|560|600|foo4|10|-
chr1|600|700|foo5|10|+
'''.replace('|', '\t')))


class StatsTest(unittest.TestCase):
    def testBedStats(self):
        stats = ngsutils.bed.stats.BedStats(bedtest)
        self.assertEqual(stats.total, 5)
        self.assertEqual(stats.size, 300)
        self.assertEqual(stats.lengths.mean(), 60.0)
        self.assertEqual(stats.refs['chr1'], 5)

    def testBedStatsGTF(self):
        ''' MISSING TEST / GTF '''
        pass

if __name__ == '__main__':
    unittest.main()
