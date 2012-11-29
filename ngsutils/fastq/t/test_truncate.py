#!/usr/bin/env python
'''
Tests for fastqutils truncate
'''

import unittest
import StringIO

import ngsutils.fastq.truncate
from ngsutils.fastq import FASTQ


class TruncateTest(unittest.TestCase):
    def testTruncate(self):
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
        ngsutils.fastq.truncate.fastq_truncate(FASTQ(fileobj=fq), 4, quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
@foo
ACGT
+
;;;;
@bar comment
ACGT
+
;;;;
''')

    def testTruncateCS(self):
        fq = StringIO.StringIO('''\
@foo
T01230123
+
;;;;;;;;
@bar comment
T01230123
+
;;;;;;;;
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.truncate.fastq_truncate(FASTQ(fileobj=fq), 4, quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
@foo
T0123
+
;;;;
@bar comment
T0123
+
;;;;
''')

if __name__ == '__main__':
    unittest.main()
