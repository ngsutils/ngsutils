#!/usr/bin/env python
'''
Tests for bamutils split
'''

import os
import unittest

import ngsutils.bam
import ngsutils.bam.merge


class MergeTest(unittest.TestCase):
    def testMerge(self):
        '''
        Merge test.bam and test2.bam to tmp.bam
        '''
        fname1 = os.path.join(os.path.dirname(__file__), 'test.bam')
        fname2 = os.path.join(os.path.dirname(__file__), 'test2.bam')
        outfname = os.path.join(os.path.dirname(__file__), 'tmp.bam')

        ngsutils.bam.merge.bam_merge(outfname, [fname1, fname2], quiet=True)

        bam = ngsutils.bam.bam_open(outfname)
        for read in ngsutils.bam.bam_iter(bam):
            if read.qname == 'A':  # AS chr2
                self.assertEqual(bam.getrname(read.tid), 'chr2')
            elif read.qname == 'B':  # AS tie
                self.assertEqual(bam.getrname(read.tid), 'chr1')
            elif read.qname == 'C':  # AS tie
                self.assertEqual(bam.getrname(read.tid), 'chr1')
            elif read.qname == 'D':  # AS tie
                self.assertEqual(bam.getrname(read.tid), 'chr1')
            elif read.qname == 'E':  # AS chr2
                self.assertEqual(bam.getrname(read.tid), 'chr2')
            elif read.qname == 'F':  # AS chr1
                self.assertEqual(bam.getrname(read.tid), 'chr1')
            elif read.qname == 'Z':  # still unmapped
                self.assertTrue(read.is_unmapped)

    def tearDown(self):
        outfname = os.path.join(os.path.dirname(__file__), 'tmp.bam')
        os.unlink(outfname)

if __name__ == '__main__':
    unittest.main()
