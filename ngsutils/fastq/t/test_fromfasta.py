#!/usr/bin/env python
'''
Tests for fastqutils fromfasta
'''

import unittest
import StringIO

from ngsutils.support import FASTA
import ngsutils.fastq.fromfasta


class FromFASTATest(unittest.TestCase):
    def testConstant(self):
        fa = StringIO.StringIO('''\
>foo_F3
ACTGACTG
''')
        fasta = FASTA(fileobj=fa)
        out = StringIO.StringIO('')
        ngsutils.fastq.fromfasta.merge_files(fasta, None, common_qual=';', out=out)

        self.assertEqual(out.getvalue(), '''\
@foo
ACTGACTG
+
;;;;;;;;
''')

    def testSuffix(self):
        fa = StringIO.StringIO('''\
>foo_F3
ACTGACTG
''')
        fasta = FASTA(fileobj=fa)
        out = StringIO.StringIO('')
        ngsutils.fastq.fromfasta.merge_files(fasta, None, common_qual=';', suffix='_bar', out=out)

        self.assertEqual(out.getvalue(), '''\
@foo_bar
ACTGACTG
+
;;;;;;;;
''')

    def testError(self):
        fa = StringIO.StringIO('''\
>foo_F3
ACTGACTG
''')
        qa = StringIO.StringIO('''\
>bar_F3
10 10 10 10 20 20 20 20
''')
        fasta = FASTA(fileobj=fa)
        qual = FASTA(fileobj=qa)
        out = StringIO.StringIO('')
        self.assertRaises(ValueError, ngsutils.fastq.fromfasta.merge_files, *(fasta, qual), **{'out': out})

    def testError2(self):
        fa = StringIO.StringIO('''\
>foo_F3
ACTGACTG
''')
        qa = StringIO.StringIO('''\
>foo_F3
10 10 10 20 20 20 20 A
''')
        fasta = FASTA(fileobj=fa)
        qual = FASTA(fileobj=qa)
        out = StringIO.StringIO('')
        self.assertRaises(ValueError, ngsutils.fastq.fromfasta.merge_files, *(fasta, qual), **{'out': out})

    def testQualFile(self):
        fa = StringIO.StringIO('''\
>foo_F3
ACTGACTG
''')
        qa = StringIO.StringIO('''\
>foo_F3
10 10 10 10 0 93 -5 99
''')
        fasta = FASTA(fileobj=fa)
        qual = FASTA(fileobj=qa, qual=True)
        out = StringIO.StringIO('')
        ngsutils.fastq.fromfasta.merge_files(fasta, qual, out=out)

        self.assertEqual(out.getvalue(), '''\
@foo
ACTGACTG
+
++++!~!~
''')

    def testQualMultiLine(self):
        fa = StringIO.StringIO('''\
>foo_F3
ACTG
ACTG
''')
        qa = StringIO.StringIO('''\
>foo_F3
10 10 10 10
0 93 -5 99
''')
        fasta = FASTA(fileobj=fa)
        qual = FASTA(fileobj=qa, qual=True)
        out = StringIO.StringIO('')
        ngsutils.fastq.fromfasta.merge_files(fasta, qual, out=out)

        self.assertEqual(out.getvalue(), '''\
@foo
ACTGACTG
+
++++!~!~
''')

    def testColor(self):
        fa = StringIO.StringIO('''\
>foo_F3
T01230123
''')
        qa = StringIO.StringIO('''\
>foo_F3
10 10 10 10 20 20 20 20
''')
        fasta = FASTA(fileobj=fa)
        qual = FASTA(fileobj=qa)
        out = StringIO.StringIO('')
        ngsutils.fastq.fromfasta.merge_files(fasta, qual, out=out)

        self.assertEqual(out.getvalue(), '''\
@foo
T01230123
+
++++5555
''')

    def testConstantColor(self):
        fa = StringIO.StringIO('''\
>foo_F3
T01230123
''')
        fasta = FASTA(fileobj=fa)
        out = StringIO.StringIO('')
        ngsutils.fastq.fromfasta.merge_files(fasta, None, common_qual=';', out=out)

        self.assertEqual(out.getvalue(), '''\
@foo
T01230123
+
;;;;;;;;
''')

if __name__ == '__main__':
    unittest.main()
