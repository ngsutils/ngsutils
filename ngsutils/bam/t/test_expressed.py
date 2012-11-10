#!/usr/bin/env python
'''
Tests for bamutils expressed
'''

import os
import StringIO
import unittest

import ngsutils.bam
import ngsutils.bam.expressed

infname = os.path.join(os.path.dirname(__file__), 'test4.bam')


class ExpressedTest(unittest.TestCase):
    def testExpressed(self):
        self._run_test(['chr1|100|259|+', 'chr1|500|570|+'])

    def testExpressedNoStrand(self):
        self._run_test(['chr1|100|259|+', 'chr1|500|570|+', 'chr1|300|425|+'], nostrand=True)

    def testExpressedSingle(self):
        self._run_test(['chr1|100|259|+', 'chr1|500|570|+', 'chr1|300|350|+',
                        'chr1|375|425|+', 'chr1|700|750|+', 'chr1|800|850|+',
                        'chr1|861|911|+', 'chr1|320|370|-'], min_read_count=1)

    def testExpressedUniq(self):
        self._run_test(['chr1|100|259|+', 'chr1|500|570|+'], only_uniq_starts=True)

    def testExpressedUniqNoStrand(self):
        self._run_test(['chr1|100|259|+', 'chr1|500|570|+', 'chr1|300|425|+'], only_uniq_starts=True, nostrand=True)

    def testExpressedLongMerge(self):
        self._run_test(['chr1|100|259|+', 'chr1|500|570|+', 'chr1|800|911|+'], merge_distance=11)

    def _run_test(self, valid, *args, **kwargs):
        out = StringIO.StringIO('')
        ngsutils.bam.expressed.bam_find_regions(infname, out=out, *args, **kwargs)
        checks = [False, ] * len(valid)
        extras = False

        for line in out.getvalue().split('\n'):
            if line:
                cols = line.strip().split('\t')
                k = '%s|%s|%s|%s' % (cols[0], cols[1], cols[2], cols[5])
                found = False
                for i, v in enumerate(valid):
                    if v == k:
                        checks[i] = True
                        found = True

                if not found:
                    extras = True
                    print "\nExtra: %s" % line.strip()

        for check in checks:
            self.assertTrue(check)
        self.assertFalse(extras)


if __name__ == '__main__':
    unittest.main()
