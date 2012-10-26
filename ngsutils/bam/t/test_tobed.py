#!/usr/bin/env python
'''
Tests for bamutils tag
'''

import unittest
from mock import *

import ngsutils.bam.tobed
import StringIO


class ToBEDTest(unittest.TestCase):
    def testBEDWriter(self):
        sio = StringIO.StringIO("")

        read = Mock()
        read.is_unmapped = False
        read.qname = 'readname'
        read.pos = 100
        read.aend = 150
        read.is_reverse = False

        ngsutils.bam.tobed.write_read(read, 'chrom', sio)
        self.assertEqual(sio.getvalue(), 'chrom\t100\t150\treadname\t0\t+\n')
        sio.close()

        read.is_reverse = True

        sio = StringIO.StringIO("")
        ngsutils.bam.tobed.write_read(read, 'chrom', sio)
        self.assertEqual(sio.getvalue(), 'chrom\t100\t150\treadname\t0\t-\n')
        sio.close()

        read.is_unmapped = True

        sio = StringIO.StringIO("")
        ngsutils.bam.tobed.write_read(read, 'chrom', sio)
        self.assertEqual(sio.getvalue(), '')
        sio.close()

if __name__ == '__main__':
    unittest.main()
