#!/usr/bin/env python
'''
Tests for bamutils split
'''

import os
import unittest

import ngsutils.bam
import ngsutils.bam.split


class SplitTest(unittest.TestCase):
    def testSplit(self):
        '''
        Merge test.bam and test2.bam to tmp.bam
        '''
        fname1 = os.path.join(os.path.dirname(__file__), 'test.bam')
        outfname = os.path.join(os.path.dirname(__file__), 'tmp')
        outfname1 = os.path.join(os.path.dirname(__file__), 'tmp.1.bam')
        outfname2 = os.path.join(os.path.dirname(__file__), 'tmp.2.bam')

        ngsutils.bam.split.bam_split(fname1, outfname, 4)

        self.assertTrue(os.path.exists(outfname1))
        self.assertTrue(os.path.exists(outfname2))

        foundA = False
        foundB = False
        foundC = False
        foundE = False
        foundOther = False

        bam = ngsutils.bam.bam_open(outfname1)
        for read in ngsutils.bam.bam_iter(bam):
            if read.qname == 'A':
                foundA = True
            elif read.qname == 'B':
                foundB = True
            elif read.qname == 'C':
                foundC = True
            elif read.qname == 'E':
                foundE = True
            else:
                foundOther = True

        self.assertTrue(foundA)
        self.assertTrue(foundB)
        self.assertTrue(foundC)
        self.assertTrue(foundE)
        self.assertFalse(foundOther)

    def tearDown(self):
        outfname1 = os.path.join(os.path.dirname(__file__), 'tmp.1.bam')
        outfname2 = os.path.join(os.path.dirname(__file__), 'tmp.2.bam')

        os.unlink(outfname1)
        os.unlink(outfname2)

if __name__ == '__main__':
    unittest.main()
