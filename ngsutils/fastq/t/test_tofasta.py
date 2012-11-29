#!/usr/bin/env python
'''
Tests for fastqutils tofasta
'''

import unittest
import StringIO

import ngsutils.fastq.tofasta
from ngsutils.fastq import FASTQ


class ToFASTATest(unittest.TestCase):
    def testSeq(self):
        fq = StringIO.StringIO('''\
@foo
ACGTACGT
+
;;;;;;;;
@bar comment
ACGTACGT
+
;;;;;;;;
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.tofasta.export_fasta(FASTQ(fileobj=fq), qual=False, quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
>foo
ACGTACGT
>bar comment
ACGTACGT
''')

    def testQual(self):
        fq = StringIO.StringIO('''\
@foo
ACGTACGT
+
!AI;;;;;
@bar comment
ACGTACGT
+
;;;;;;;;
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.tofasta.export_fasta(FASTQ(fileobj=fq), qual=True, quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
>foo
0 32 40 26 26 26 26 26
>bar comment
26 26 26 26 26 26 26 26
''')


if __name__ == '__main__':
    unittest.main()
