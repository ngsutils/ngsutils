#!/usr/bin/env python
'''
Tests for bamutils export
'''

import StringIO
import unittest

import ngsutils.bam
import ngsutils.bam.export

from ngsutils.bam.t import MockBam


testbam1 = MockBam(['chr1'])
testbam1.add_read('foo1', 'atcgatcg', 'AAAAAAAA', 0, 0, 8, '8M')
testbam1.add_read('foo2', 'atcgtttt', 'AAAABBBB', 0, 4, 12, '8M', is_reverse=True)
testbam1.add_read('foo3', 'gggggggg', 'CCCCCCCC')

testbam2 = MockBam(['chr1'])
testbam2.add_read('bar1', 'atcgatcg', 'AAAAAAAA', 0, 0, 8, '8M', rnext=0, pnext=100, is_paired=True, is_read1=True)
testbam2.add_read('bar1', 'gctagcta', 'AAAAAAAA', 0, 100, 108, '8M', is_paired=True, is_read2=True)
testbam2.add_read('bar2', 'atcgatcg', 'AAAAAAAA', is_paired=True, is_read1=True)
testbam2.add_read('bar2', 'gctagcta', 'AAAAAAAA', is_paired=True, is_read2=True)

testbam3 = MockBam(['chr1'])
testbam3.add_read('baz1', 'atcgatcg', 'AAAAAAAA', 0, 0, 8, '8M', mapq=10, isize=100, tlen=120, tags=[('ZZ', 'foo'), ('ZY', 100), ('ZX', 1.0), ('ZW', 'a'), ])


class ExportTest(unittest.TestCase):
    def testExport1(self):
        self._run_test(testbam1, ['foo1|chr1|1|8M', 'foo2|chr1|5|8M', 'foo3|*|0|*'], fields=['-name', '-ref', '-pos', '-cigar'])

    def testExportMapped(self):
        self._run_test(testbam1, ['foo1|chr1|1|8M', 'foo2|chr1|5|8M'], fields=['-name', '-ref', '-pos', '-cigar'], unmapped=False)

    def testExportUnmapped(self):
        self._run_test(testbam1, ['foo3|*|0|*'], fields=['-name', '-ref', '-pos', '-cigar'], mapped=False)

    def testExportWhite(self):
        self._run_test(testbam1, ['foo2|chr1|5|8M', ], fields=['-name', '-ref', '-pos', '-cigar'], whitelist=['foo2'])

    def testExportBlack(self):
        self._run_test(testbam1, ['foo1|chr1|1|8M', 'foo3|*|0|*'], fields=['-name', '-ref', '-pos', '-cigar'], blacklist=['foo2'])

    def testSeqQual(self):
        self._run_test(testbam1, ['ATCGATCG|AAAAAAAA', 'ATCGTTTT|AAAABBBB', 'GGGGGGGG|CCCCCCCC'], fields=['-seq', '-qual'])

    def testFlags(self):
        self._run_test(testbam1, ['foo1|0', 'foo2|16', 'foo3|4'], fields=['-name', '-flags'])
        self._run_test(testbam2, ['bar1|65', 'bar1|137', 'bar2|77', 'bar2|141'], fields=['-name', '-flags'])

    def testPairs(self):
        self._run_test(testbam2, ['bar1|65|chr1|chr1|1|101', 'bar1|137|chr1|*|101|0'], fields=['-name', '-flags', '-ref', '-nextref', '-pos', '-nextpos'], unmapped=False)
        self._run_test(testbam2, ['bar2|77|*|*|0|0', 'bar2|141|*|*|0|0'], fields=['-name', '-flags', '-ref', '-nextref', '-pos', '-nextpos'], mapped=False)

    def testExtended(self):
        self._run_test(testbam3, ['baz1|10|100|120'], fields=['-name', '-mapq', '-isize', '-tlen'])

    def testTags(self):
        self._run_test(testbam3, ['baz1|ZZ:Z:foo|ZY:i:100|ZX:f:1.0|ZW:A:a'], fields=['-name', '-tag:*'])
        self._run_test(testbam3, ['baz1|foo'], fields=['-name', '-tag:ZZ'])


    def _run_test(self, testbam, valid, *args, **kwargs):
        out = StringIO.StringIO('')
        ngsutils.bam.export.bam_export(testbam, out=out, quiet=True, *args, **kwargs)
        checks = [False, ] * len(valid)
        extras = False

        for line in out.getvalue().split('\n'):
            if line:
                cols = line.strip().split('\t')
                k = '|'.join([str(x) for x in cols])
                found = False
                for i, v in enumerate(valid):
                    if v == k:
                        checks[i] = True
                        found = True

                if not found:
                    extras = True
                    print "Extra: %s" % line.strip()

        for check in checks:
            self.assertTrue(check)
        self.assertFalse(extras)


if __name__ == '__main__':
    unittest.main()
