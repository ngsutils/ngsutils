#!/usr/bin/env python
'''
Tests for fastqutils filter
'''

import unittest
import StringIO

import ngsutils.fastq.properpairs
from ngsutils.fastq import FASTQ


class ProperPairsTest(unittest.TestCase):
    def testMissing1(self):
        fq1 = StringIO.StringIO('''\
@foo comment
ACGTACGT
+
;;;;;;;;
@bar
ACGTACGT
+
;;;;;;;;
@baz
ACGTACGT
+
;;;;;;;;
@quux
ACGTACGT
+
;;;;;;;;
''')
        fq2 = StringIO.StringIO('''\
@foo comment
ACGTACGT
+
;;;;;;;;
@bar
ACGTACGT
+
;;;;;;;;
@quux
ACGTACGT
+
;;;;;;;;
''')

        out1 = StringIO.StringIO('')
        out2 = StringIO.StringIO('')

        fastq1 = FASTQ(fileobj=fq1)
        fastq2 = FASTQ(fileobj=fq2)

        ngsutils.fastq.properpairs.find_fastq_pairs(fastq1, fastq2, out1, out2, quiet=True)
        self.assertEqual(out1.getvalue(), out2.getvalue())
        self.assertEqual(out1.getvalue(), fq2.getvalue())

    def testMissing2(self):
        fq1 = StringIO.StringIO('''\
@foo comment
ACGTACGT
+
;;;;;;;;
@bar
ACGTACGT
+
;;;;;;;;
@quux
ACGTACGT
+
;;;;;;;;
''')
        fq2 = StringIO.StringIO('''\
@foo comment
ACGTACGT
+
;;;;;;;;
@bar
ACGTACGT
+
;;;;;;;;
@baz
ACGTACGT
+
;;;;;;;;
@quux
ACGTACGT
+
;;;;;;;;
''')

        out1 = StringIO.StringIO('')
        out2 = StringIO.StringIO('')

        fastq1 = FASTQ(fileobj=fq1)
        fastq2 = FASTQ(fileobj=fq2)

        ngsutils.fastq.properpairs.find_fastq_pairs(fastq1, fastq2, out1, out2, quiet=True)
        self.assertEqual(out1.getvalue(), out2.getvalue())
        self.assertEqual(out1.getvalue(), fq1.getvalue())

    def testMissingBoth(self):
        fq1 = StringIO.StringIO('''\
@foo comment
ACGTACGT
+
;;;;;;;;
@bar
ACGTACGT
+
;;;;;;;;
@quux
ACGTACGT
+
;;;;;;;;
''')
        fq2 = StringIO.StringIO('''\
@foo comment
ACGTACGT
+
;;;;;;;;
@baz
ACGTACGT
+
;;;;;;;;
@quux
ACGTACGT
+
;;;;;;;;
''')

        out1 = StringIO.StringIO('')
        out2 = StringIO.StringIO('')

        fastq1 = FASTQ(fileobj=fq1)
        fastq2 = FASTQ(fileobj=fq2)

        ngsutils.fastq.properpairs.find_fastq_pairs(fastq1, fastq2, out1, out2, quiet=True)
        self.assertEqual(out1.getvalue(), out2.getvalue())
        self.assertEqual(out1.getvalue(), '''\
@foo comment
ACGTACGT
+
;;;;;;;;
@quux
ACGTACGT
+
;;;;;;;;
''')

    def testMissingTail(self):
        fq1 = StringIO.StringIO('''\
@foo comment
ACGTACGT
+
;;;;;;;;
@bar
ACGTACGT
+
;;;;;;;;
@baz
ACGTACGT
+
;;;;;;;;
@quux
ACGTACGT
+
;;;;;;;;
''')
        fq2 = StringIO.StringIO('''\
@foo comment
ACGTACGT
+
;;;;;;;;
@baz
ACGTACGT
+
;;;;;;;;
''')

        out1 = StringIO.StringIO('')
        out2 = StringIO.StringIO('')

        fastq1 = FASTQ(fileobj=fq1)
        fastq2 = FASTQ(fileobj=fq2)

        ngsutils.fastq.properpairs.find_fastq_pairs(fastq1, fastq2, out1, out2, quiet=True)
        self.assertEqual(out1.getvalue(), out2.getvalue())
        self.assertEqual(out1.getvalue(), fq2.getvalue())

    def testMissingHead(self):
        fq1 = StringIO.StringIO('''\
@foo comment
ACGTACGT
+
;;;;;;;;
@bar
ACGTACGT
+
;;;;;;;;
@baz
ACGTACGT
+
;;;;;;;;
@quux
ACGTACGT
+
;;;;;;;;
''')
        fq2 = StringIO.StringIO('''\
@bar
ACGTACGT
+
;;;;;;;;
@quux
ACGTACGT
+
;;;;;;;;
''')

        out1 = StringIO.StringIO('')
        out2 = StringIO.StringIO('')

        fastq1 = FASTQ(fileobj=fq1)
        fastq2 = FASTQ(fileobj=fq2)

        ngsutils.fastq.properpairs.find_fastq_pairs(fastq1, fastq2, out1, out2, quiet=True)
        self.assertEqual(out1.getvalue(), out2.getvalue())
        self.assertEqual(out1.getvalue(), fq2.getvalue())


if __name__ == '__main__':
    unittest.main()
