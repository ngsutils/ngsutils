#!/usr/bin/env python
'''
Tests for bamutils extract
'''

import unittest
import os

import ngsutils.bam
import ngsutils.bam.extract

from ngsutils.bam.t import MockBam, _matches

testbam1 = MockBam(['chr1'])
testbam1.add_read('foo1', tid=0, pos=100, aend=150, cigar='50M')  # too early
testbam1.add_read('foo2', tid=0, pos=200, aend=250, cigar='50M')  # in range, + strand
testbam1.add_read('foo3', tid=0, pos=1000, aend=1050, cigar='50M')  # too late
testbam1.add_read('foo4', tid=0, pos=200, aend=350, cigar='25M100N25M')  # in range, start, + strand
testbam1.add_read('foo5', tid=0, pos=200, aend=250, cigar='50M', is_reverse=True)  # in range, start, - strand
testbam1.add_read('foo6', tid=0, pos=2000, aend=3050, cigar='20M500N20M500M10M')  # touches (2000,2020), (2520, 2540), (3040, 3050)


class ExtractTest(unittest.TestCase):
    def setUp(self):
        self.fname1 = os.path.join(os.path.dirname(__file__), 'testbam1')
        with open(self.fname1, 'w') as f:
            f.write('chr1\t200\t250\tfoo\t1\t+\n')
            f.write('chr1\t125\t170\tfoo\t1\t+\n')

        self.fname2 = os.path.join(os.path.dirname(__file__), 'testbam2')
        with open(self.fname2, 'w') as f:
            f.write('chr1\t200\t250\n')
            f.write('chr1\t125\t170\n')

    def tearDown(self):
        os.unlink(self.fname1)
        os.unlink(self.fname2)

    def testExtractBED(self):
        outbam = MockBam(['chr1'])
        ngsutils.bam.extract.bam_extract(testbam1, outbam, self.fname1, quiet=True)
        passed = [x.qname for x in outbam]
        self.assertTrue(_matches(['foo2', 'foo1', 'foo4'], passed))

    def testExtractBED3(self):
        outbam = MockBam(['chr1'])
        ngsutils.bam.extract.bam_extract(testbam1, outbam, self.fname2, quiet=True)
        passed = [x.qname for x in outbam]
        self.assertTrue(_matches(['foo2', 'foo5', 'foo1', 'foo4'], passed))

    def testExtract(self):
        passed = [x.qname for x in ngsutils.bam.extract.bam_extract_reads(testbam1, 'chr1', 200, 250, '+')]
        self.assertTrue(_matches(['foo2', 'foo4'], passed))

        passed = [x.qname for x in ngsutils.bam.extract.bam_extract_reads(testbam1, 'chr1', 200, 250, '-')]
        self.assertTrue(_matches(['foo5'], passed))

        passed = [x.qname for x in ngsutils.bam.extract.bam_extract_reads(testbam1, 'chr1', 200, 250, None)]
        self.assertTrue(_matches(['foo2', 'foo4', 'foo5'], passed))

    def testExtractTouching(self):
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
