#!/usr/bin/env python
'''
Tests for fastqutils tag
'''

import unittest
import StringIO

import ngsutils.fastq.tag
from ngsutils.fastq import FASTQ


class TagTest(unittest.TestCase):
    def testPrefix(self):
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
        ngsutils.fastq.tag.fastq_tag(FASTQ(fileobj=fq), prefix='pre_', quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
@pre_foo
ACGTACGT
+
;;;;;;;;
@pre_bar comment
ACGTACGT
+
;;;;;;;;
''')

    def testSuffix(self):
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
        ngsutils.fastq.tag.fastq_tag(FASTQ(fileobj=fq), suffix='_post', quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
@foo_post
ACGTACGT
+
;;;;;;;;
@bar_post comment
ACGTACGT
+
;;;;;;;;
''')

    def testPrefixSuffix(self):
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
        ngsutils.fastq.tag.fastq_tag(FASTQ(fileobj=fq), prefix='pre_', suffix='_post', quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
@pre_foo_post
ACGTACGT
+
;;;;;;;;
@pre_bar_post comment
ACGTACGT
+
;;;;;;;;
''')

    def testAssert(self):
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
        fastq = FASTQ(fileobj=fq)
        self.assertRaises(ValueError, ngsutils.fastq.tag.fastq_tag, *[fastq], **{'quiet': True, 'out': out})

if __name__ == '__main__':
    unittest.main()
