#!/usr/bin/env python
'''
Tests for bamutils extract
'''

import unittest
import os

import ngsutils.bam
from ngsutils.bam.innerdist import bam_innerdist

from ngsutils.bam.t import MockBam, _matches

testbam1 = MockBam(['chr1'], insert_order=True)
testbam1.add_read('foo1', tid=0, pos=100, aend=150, cigar='50M')
testbam1.add_read('foo2', tid=0, pos=200, aend=250, cigar='50M')
testbam1.add_read('foo3', tid=0, pos=300, aend=350, cigar='50M')
testbam1.add_read('foo4', tid=0, pos=400, aend=450, cigar='50M')
testbam1.add_read('foo5', tid=0, pos=400, aend=450, cigar='50M')
testbam1.add_read('foo6', tid=-1)

testbam2 = MockBam(['chr1'], insert_order=True)
testbam2.add_read('foo1', tid=0, pos=200, aend=250, cigar='50M')
testbam2.add_read('foo2', tid=0, pos=300, aend=350, cigar='50M')
testbam2.add_read('foo3', tid=0, pos=450, aend=500, cigar='50M')
testbam2.add_read('foo4', tid=0, pos=550, aend=600, cigar='50M')
testbam2.add_read('foo5', tid=-1)
testbam2.add_read('foo6', tid=0, pos=500, aend=550, cigar='50M')

testbam3 = MockBam(['chr1'], insert_order=True)
testbam3.add_read('foo1', tid=0, pos=100, aend=150, cigar='50M')
testbam3.add_read('foo2', tid=0, pos=200, aend=250, cigar='50M')
testbam3.add_read('foo3', tid=0, pos=300, aend=350, cigar='50M')
testbam3.add_read('foo4', tid=0, pos=400, aend=450, cigar='50M')

testbam4 = MockBam(['chr1'], insert_order=True)
testbam4.add_read('foo1', tid=0, pos=200, aend=250, cigar='50M')
testbam4.add_read('foo3', tid=0, pos=450, aend=500, cigar='50M')
testbam4.add_read('foo2', tid=0, pos=300, aend=350, cigar='50M')
testbam4.add_read('foo4', tid=0, pos=550, aend=600, cigar='50M')

class InnerDistTest(unittest.TestCase):
    def testError(self):
        failed = False
        try:
            bam_innerdist(testbam3, testbam4)
        except ValueError:
            failed = True

        self.assertTrue(failed)


    def testDist(self):
        total, proper, mean, stdev, o_count = bam_innerdist(testbam1, testbam2)

        self.assertEqual(total, 6)
        self.assertEqual(proper, 4)
        self.assertEqual(mean, 75.0)
        self.assertEqual(round(stdev, 5), 28.86751)


if __name__ == '__main__':
    unittest.main()
