#!/usr/bin/env python
'''
Tests for bamutils tobed
'''

import unittest
from ngsutils.bam.t import MockRead

import ngsutils.bam.tobed
import StringIO


class ToBEDTest(unittest.TestCase):
    def testBEDWriter(self):
        sio = StringIO.StringIO("")

        read = MockRead('readname', tid=1, pos=100, aend=150)

        # write out correct strand
        ngsutils.bam.tobed.write_read(read, 'chrom', sio)
        self.assertEqual(sio.getvalue(), 'chrom\t100\t150\treadname\t0\t+\n')
        sio.close()

        # write out correct strand
        read.is_reverse = True
        sio = StringIO.StringIO("")
        ngsutils.bam.tobed.write_read(read, 'chrom', sio)
        self.assertEqual(sio.getvalue(), 'chrom\t100\t150\treadname\t0\t-\n')
        sio.close()

        # skip unmapped reads
        read.is_unmapped = True
        sio = StringIO.StringIO("")
        ngsutils.bam.tobed.write_read(read, 'chrom', sio)
        self.assertEqual(sio.getvalue(), '')
        sio.close()

if __name__ == '__main__':
    unittest.main()
