#!/usr/bin/env python
'''
Tests for fastqutils names
'''

import unittest
import StringIO

import ngsutils.fastq.names
from ngsutils.fastq import FASTQ


class NamesTest(unittest.TestCase):
    def testNames(self):
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
        ngsutils.fastq.names.export_names(FASTQ(fileobj=fq), out=out, quiet=True)
        self.assertEqual(out.getvalue(), 'foo\nbar\n')

if __name__ == '__main__':
    unittest.main()
