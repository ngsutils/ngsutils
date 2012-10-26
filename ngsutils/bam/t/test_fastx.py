#!/usr/bin/env python
'''
Tests for bamutils tag
'''

import unittest
from mock import *

import ngsutils.bam.tofastq
import StringIO

read = Mock()
read.qname = 'foo'
read.seq = 'ACGT'
read.qual = 'AAAA'
read.tags = [('CS', 'T0123'), ('CQ', 'BBBB')]


def opt(tag):
    for k, v in read.tags:
        if k == tag:
            return v
    return None
read.opt = Mock(side_effect=opt)


class FASTXTest(unittest.TestCase):
    def testFASTQ(self):
        sio = StringIO.StringIO("")
        ngsutils.bam.tofastq.write_fastq(read, sio)

        self.assertEqual(sio.getvalue(), '''\
@foo
ACGT
+
AAAA
''')

    def testFASTA(self):
        sio = StringIO.StringIO("")
        ngsutils.bam.tofastq.write_fasta(read, sio)

        self.assertEqual(sio.getvalue(), '''\
>foo
ACGT
''')

    def testFASTQ_cs(self):
        sio = StringIO.StringIO("")
        ngsutils.bam.tofastq.write_fastq(read, sio, colorspace=True)

        self.assertEqual(sio.getvalue(), '''\
@foo
T0123
+
BBBB
''')

    def testFASTA_cs(self):
        sio = StringIO.StringIO("")
        ngsutils.bam.tofastq.write_fasta(read, sio, colorspace=True)

        self.assertEqual(sio.getvalue(), '''\
>foo
T0123
''')


if __name__ == '__main__':
    unittest.main()
