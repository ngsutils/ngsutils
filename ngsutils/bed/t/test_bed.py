#!/usr/bin/env python
'''
Tests for bedutils
'''

import unittest
import os

from ngsutils.bam.t import _matches
from ngsutils.bed import BedFile


class BedTest(unittest.TestCase):
    def testBed(self):
        fname = os.path.join(os.path.dirname(__file__), 'test.bed')
        valid = ['chr1|100|150|foo|1|+',
                   'chr1|100|150|foo|1|-',
                   'chr1|200|250|foo|1|+',
                   'chr1|300|350|foo|1|-',
                   ]

        regions = ['%s|%s|%s|%s|%s|%s' % (x.chrom, x.start, x.end, x.name, x.score, x.strand) for x in BedFile(fname)]
        self.assertTrue(_matches(valid, regions))


if __name__ == '__main__':
    unittest.main()
