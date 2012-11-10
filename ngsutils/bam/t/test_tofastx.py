#!/usr/bin/env python
'''
Tests for bamutils tofasta / tofastq
'''

import unittest
from ngsutils.bam.t import MockRead

import ngsutils.bam.tofastq
import StringIO

read = MockRead('foo', 'ACGT', 'AAAA', tags=[('CS', 'T0123'), ('CQ', 'BBBB')])


class FASTXTest(unittest.TestCase):
    def testFASTQ(self):
        sio = StringIO.StringIO("")
        ngsutils.bam.tofastq.write_fastq(read, sio)
        self.assertEqual(sio.getvalue(), '@foo\nACGT\n+\nAAAA\n')
        sio.close()

    def testFASTA(self):
        sio = StringIO.StringIO("")
        ngsutils.bam.tofastq.write_fasta(read, sio)
        self.assertEqual(sio.getvalue(), '>foo\nACGT\n')
        sio.close()

    def testFASTQ_cs(self):
        sio = StringIO.StringIO("")
        ngsutils.bam.tofastq.write_fastq(read, sio, colorspace=True)
        self.assertEqual(sio.getvalue(), '@foo\nT0123\n+\nBBBB\n')
        sio.close()

    def testFASTA_cs(self):
        sio = StringIO.StringIO("")
        ngsutils.bam.tofastq.write_fasta(read, sio, colorspace=True)
        self.assertEqual(sio.getvalue(), '>foo\nT0123\n')
        sio.close()


if __name__ == '__main__':
    unittest.main()
