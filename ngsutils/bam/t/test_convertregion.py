#!/usr/bin/env python
'''
Tests for bamutils convertregion
'''

import os
import unittest

import ngsutils.bam
import ngsutils.bam.convertregion

infname = os.path.join(os.path.dirname(__file__), 'test3.bam')
outfname = os.path.join(os.path.dirname(__file__), 'tmp-convert.bam')
chromsizes = os.path.join(os.path.dirname(__file__), 'test.chromsizes.txt')


class ConvertRegionTest(unittest.TestCase):
    def testConvertRegion(self):
        '''
        Convert region/junction coordinates
        '''
        ngsutils.bam.convertregion.bam_convertregion(infname, outfname, chromsizes, quiet=True)

        bam = ngsutils.bam.bam_open(outfname)
        for read in ngsutils.bam.bam_iter(bam):
            if read.qname == 'A':
                self.assertEqual(bam.getrname(read.tid), 'chr1')
                self.assertEqual(read.pos, 109)  # 0-based in BAM, 1-based in SAM
                self.assertEqual(read.cigar, ngsutils.bam.cigar_fromstr('50M'))
            elif read.qname == 'B':
                self.assertEqual(bam.getrname(read.tid), 'chr1')
                self.assertEqual(read.pos, 174)
                self.assertEqual(read.cigar, ngsutils.bam.cigar_fromstr('26M100N24M'))
            elif read.qname == 'C':
                self.assertEqual(bam.getrname(read.tid), 'chr1')
                self.assertEqual(read.pos, 197)
                self.assertEqual(read.cigar, ngsutils.bam.cigar_fromstr('3M100N47M'))
            elif read.qname == 'Z':
                self.assertTrue(read.is_unmapped)

    def testConvertRegionOverlap(self):
        '''
        Convert region/junction coordinates, require overlap with junction
        '''
        ngsutils.bam.convertregion.bam_convertregion(infname, outfname, chromsizes, enforce_overlap=True, quiet=True)
        foundA = False
        foundB = False
        foundC = False
        foundZ = False

        bam = ngsutils.bam.bam_open(outfname)
        for read in ngsutils.bam.bam_iter(bam):
            if read.qname == 'A':
                foundA = True
            elif read.qname == 'B':
                self.assertEqual(bam.getrname(read.tid), 'chr1')
                self.assertEqual(read.pos, 174)
                self.assertEqual(read.cigar, ngsutils.bam.cigar_fromstr('26M100N24M'))
                foundB = True
            elif read.qname == 'C':
                foundC = True
            elif read.qname == 'Z':
                foundZ = True

        self.assertFalse(foundA)
        self.assertTrue(foundB)
        self.assertFalse(foundC)
        self.assertFalse(foundZ)

    def tearDown(self):
        outfname = os.path.join(os.path.dirname(__file__), 'tmp-convert.bam')
        os.unlink(outfname)

if __name__ == '__main__':
    unittest.main()
