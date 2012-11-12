#!/usr/bin/env python
'''
Tests for bamutils extract
'''

import unittest

import ngsutils.bam
import ngsutils.bam.extract

from ngsutils.bam.t import MockBam

testbam1 = MockBam(['chr1'])
testbam1.add_read('foo1', tid=0, pos=100, aend=150, cigar='50M')  # too early
testbam1.add_read('foo2', tid=0, pos=200, aend=250, cigar='50M')  # in range, + strand
testbam1.add_read('foo3', tid=0, pos=1000, aend=1050, cigar='50M')  # too late
testbam1.add_read('foo4', tid=0, pos=200, aend=350, cigar='25M100N25M')  # in range, start, + strand
testbam1.add_read('foo5', tid=0, pos=200, aend=250, cigar='50M', is_reverse=True)  # in range, start, - strand
testbam1.add_read('foo6', tid=0, pos=2000, aend=3050, cigar='20M500N20M500M10M')  # touches (2000,2020), (2520, 2540), (3040, 3050)


def _matches(valid, queries):
    extra = False
    check = [False, ] * len(valid)
    for q in queries:
        found = False
        for i, v in enumerate(valid):
            if v == q:
                check[i] = True
                found = True
                break
        if not found:
            extra = True

    if not check or extra:
        return False

    return True


class ExtractTest(unittest.TestCase):
    def testExtract(self):
        passed = [x.qname for x in ngsutils.bam.extract.bam_extract_reads(testbam1, 'chr1', 200, 250, '+')]
        self.assertTrue(_matches(['foo2', 'foo4'], passed))

        passed = [x.qname for x in ngsutils.bam.extract.bam_extract_reads(testbam1, 'chr1', 200, 250, '-')]
        self.assertTrue(_matches(['foo5'], passed))

        passed = [x.qname for x in ngsutils.bam.extract.bam_extract_reads(testbam1, 'chr1', 200, 250, None)]
        self.assertTrue(_matches(['foo2', 'foo4', 'foo5'], passed))

    def testTouching(self):
        passed = [x.qname for x in ngsutils.bam.extract.bam_extract_reads(testbam1, 'chr1', 2520, 2540, None)]
        self.assertTrue(_matches(['foo6'], passed))

        passed = [x.qname for x in ngsutils.bam.extract.bam_extract_reads(testbam1, 'chr1', 2510, 2530, None)]
        self.assertTrue(_matches(['foo6'], passed))

        passed = [x.qname for x in ngsutils.bam.extract.bam_extract_reads(testbam1, 'chr1', 2530, 2550, None)]
        self.assertTrue(_matches(['foo6'], passed))

        passed = [x.qname for x in ngsutils.bam.extract.bam_extract_reads(testbam1, 'chr1', 2500, 2600, None)]
        self.assertTrue(_matches(['foo6'], passed))

        passed = [x.qname for x in ngsutils.bam.extract.bam_extract_reads(testbam1, 'chr1', 2530, 2535, None)]
        self.assertTrue(_matches(['foo6'], passed))

if __name__ == '__main__':
    unittest.main()
