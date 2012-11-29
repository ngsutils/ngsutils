#!/usr/bin/env python
'''
Tests for fastqutils merge
'''

import unittest
import StringIO

import ngsutils.fastq.merge
from ngsutils.fastq import FASTQ


class FASTQMergeTest(unittest.TestCase):
    def testMerge(self):
        fq1 = StringIO.StringIO('''\
@foo
ACGTACGT
+
;;;;;;;;
@bar comment
ACGTACGT
+
;;;;;;;;
''')
        fq2 = StringIO.StringIO('''\
@foo
acgtacgt
+
AAAAAAAA
@bar comment
acgtacgt
+
AAAAAAAA
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.merge.fastq_merge([FASTQ(fileobj=fq1), FASTQ(fileobj=fq2)], out=out, quiet=True)
        self.assertEqual(out.getvalue(), '''\
@foo
ACGTACGT
+
;;;;;;;;
@foo
acgtacgt
+
AAAAAAAA
@bar comment
ACGTACGT
+
;;;;;;;;
@bar comment
acgtacgt
+
AAAAAAAA
''')

    def testSplit(self):
        fq1 = StringIO.StringIO('''\
@foo/1 comment1
ACGTACGT
+
;;;;;;;;
''')
        fq2 = StringIO.StringIO('''\
@foo/2 comment2
acgtacgt
+
AAAAAAAA
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.merge.fastq_merge([FASTQ(fileobj=fq1), FASTQ(fileobj=fq2)], split_slashes=True, out=out, quiet=True)
        self.assertEqual(out.getvalue(), '''\
@foo /1 comment1
ACGTACGT
+
;;;;;;;;
@foo /2 comment2
acgtacgt
+
AAAAAAAA
''')

    def testAssert(self):
        fq1 = StringIO.StringIO('''\
@foo
ACGTACGT
+
;;;;;;;;
''')
        fq2 = StringIO.StringIO('''\
@bar
acgtacgt
+
AAAAAAAA
''')
        out = StringIO.StringIO('')
        self.assertRaises(ValueError, ngsutils.fastq.merge.fastq_merge, *[[FASTQ(fileobj=fq1), FASTQ(fileobj=fq2)], ], **{'out': out, 'quiet': True})


if __name__ == '__main__':
    unittest.main()
