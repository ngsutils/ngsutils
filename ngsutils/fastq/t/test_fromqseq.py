#!/usr/bin/env python
'''
Tests for fastqutils fromqseq
'''

import unittest
import StringIO

import ngsutils.fastq.fromqseq


class FromQseqTest(unittest.TestCase):
    def testQseq(self):
        qseq = StringIO.StringIO('''\
foo|1|2|3|4|5|0|1|ACGTACGT|hhhhhhhh|1
'''.replace('|', '\t'))
        out = StringIO.StringIO('')
        ngsutils.fastq.fromqseq.read_illumina_export(ngsutils.fastq.fromqseq.qseq_reader(fileobj=qseq), out=out)
        self.assertEqual(out.getvalue(), '''\
@foo:1:2:3:4:5 #0/1
ACGTACGT
+
IIIIIIII
''')

    def testQseqTag(self):
        qseq = StringIO.StringIO('''\
foo|1|2|3|4|5|0|1|ACGTACGT|hhhhhhhh|1
foo|1|2|3|4|6|0|1|ACGTACGT|hhhhhhhh|0
'''.replace('|', '\t'))
        out = StringIO.StringIO('')
        ngsutils.fastq.fromqseq.read_illumina_export(ngsutils.fastq.fromqseq.qseq_reader(fileobj=qseq), tag='bar_', out=out)
        self.assertEqual(out.getvalue(), '''\
@bar_foo:1:2:3:4:5 #0/1
ACGTACGT
+
IIIIIIII
''')

    def testQseqLength(self):
        qseq = StringIO.StringIO('''\
foo|1|2|3|4|5|0|1|ACGTACGT|hhhhhhhh|1
foo|1|2|3|4|6|0|1|ACGTACGT|hhhhhhhB|1
foo|1|2|3|4|7|0|1|ACGTACGT|BBBBBBBB|0
'''.replace('|', '\t'))
        out = StringIO.StringIO('')
        results = ngsutils.fastq.fromqseq.read_illumina_export(ngsutils.fastq.fromqseq.qseq_reader(fileobj=qseq), trim=True, min_length=8, out=out)
        self.assertEqual(results.reads, 3)
        self.assertEqual(results.qcfailed, 1)
        self.assertEqual(results.lenfailed, 1)
        self.assertEqual(results.passed, 1)
        self.assertEqual(results.lengths[8], 1)
        self.assertEqual(out.getvalue(), '''\
@foo:1:2:3:4:5 #0/1
ACGTACGT
+
IIIIIIII
''')

    def testQseqQCFail(self):
        qseq = StringIO.StringIO('''\
foo|1|2|3|4|5|0|1|ACGTACGT|hhhhhhhh|1
foo|1|2|3|4|6|0|1|ACGTACGT|hhhhhhhh|0
'''.replace('|', '\t'))
        out = StringIO.StringIO('')
        ngsutils.fastq.fromqseq.read_illumina_export(ngsutils.fastq.fromqseq.qseq_reader(fileobj=qseq), out=out)
        self.assertEqual(out.getvalue(), '''\
@foo:1:2:3:4:5 #0/1
ACGTACGT
+
IIIIIIII
''')

    def testQseqQC(self):
        qseq = StringIO.StringIO('''\
foo|1|2|3|4|5|0|1|ACGTACGT|hhhhhhhh|1
foo|1|2|3|4|6|0|1|ACGTACGT|hhhhhhhh|0
'''.replace('|', '\t'))
        out = StringIO.StringIO('')
        ngsutils.fastq.fromqseq.read_illumina_export(ngsutils.fastq.fromqseq.qseq_reader(fileobj=qseq), qc_remove=False, out=out)
        self.assertEqual(out.getvalue(), '''\
@foo:1:2:3:4:5 #0/1
ACGTACGT
+
IIIIIIII
@foo:1:2:3:4:6 #0/1 QCFAIL
ACGTACGT
+
IIIIIIII
''')

    def testQseqSolexa(self):
        qseq = StringIO.StringIO('''\
foo|1|2|3|4|5|0|1|ACGTACGT|@@@@hhhh|1
'''.replace('|', '\t'))
        out = StringIO.StringIO('')
        ngsutils.fastq.fromqseq.read_illumina_export(ngsutils.fastq.fromqseq.qseq_reader(fileobj=qseq), solexa_quals=True, out=out)

        # at high quals, the two scales converge (h/I -> 40/40,
        # but at lower scales (@/$ -> 0/3) they diverge
        self.assertEqual(out.getvalue(), '''\
@foo:1:2:3:4:5 #0/1
ACGTACGT
+
$$$$IIII
''')


if __name__ == '__main__':
    unittest.main()
