#!/usr/bin/env python
'''
Tests for fastqutils / docutils
'''

import os
import unittest
import doctest
import StringIO

import ngsutils.fastq
import ngsutils.fastq.fromfasta


class FASTQTest(unittest.TestCase):
    def testSimpleFile(self):
                fastq = ngsutils.fastq.FASTQ(os.path.join(os.path.dirname(__file__), 'test.fastq'))
                names = [x.name.split()[0] for x in fastq.fetch(quiet=True)]
                self.assertEqual(names, ['foo', 'foo', 'bar', 'bar', 'baz', 'baz'])

    def testSimple(self):
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

        fastq = ngsutils.fastq.FASTQ(fileobj=fq)
        self.assertEqual(fastq.is_paired, False)
        self.assertEqual(fastq.is_colorspace, False)

    def testPaired(self):
        fq = StringIO.StringIO('''\
@foo /1
ACGTACGT
+
;;;;;;;;
@foo /2 comment
ACGTACGT
+
;;;;;;;;
''')

        fastq = ngsutils.fastq.FASTQ(fileobj=fq)
        self.assertEqual(fastq.is_paired, True)
        self.assertEqual(fastq.is_colorspace, False)

    def testColorspace(self):
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

        fastq = ngsutils.fastq.FASTQ(fileobj=fq)
        self.assertEqual(fastq.is_paired, False)
        self.assertEqual(fastq.is_colorspace, True)


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(ngsutils.fastq))
    tests.addTests(doctest.DocTestSuite(ngsutils.fastq.fromfasta))
    return tests

if __name__ == '__main__':
    unittest.main()
