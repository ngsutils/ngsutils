#!/usr/bin/env python
'''
Tests for fastqutils trim
'''

import unittest
import StringIO

import ngsutils.fastq.trim
from ngsutils.fastq import FASTQ


class TrimTest(unittest.TestCase):
    def testTrim5(self):
        fq = StringIO.StringIO('''\
@foo
ACGTaaccggttccttggaa
+
;;;;;;;;;;;;;;;;;;;;
@bar - note offset
aACGTtgtgatagctacgact
+
;;;;;;;;;;;;;;;;;;;;;
@baz
AGCTtgtagatgatagataga
+
;;;;;;;;;;;;;;;;;;;;;
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.trim.fastq_trim(FASTQ(fileobj=fq), linker_5='ACGT', min_len=10, quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
@foo
aaccggttccttggaa
+
;;;;;;;;;;;;;;;;
@bar - note offset
tgtgatagctacgact
+
;;;;;;;;;;;;;;;;
@baz
AGCTtgtagatgatagataga
+
;;;;;;;;;;;;;;;;;;;;;
''')

    def testTrim3(self):
        fq = StringIO.StringIO('''\
@foo
aaccggttccttggaaACGT
+
;;;;;;;;;;;;;;;;;;;;
@bar - note offset
tgtgatagctacgactACGTa
+
;;;;;;;;;;;;;;;;;;;;;
@baz
tgtagatgatagatagaAGCT
+
;;;;;;;;;;;;;;;;;;;;;
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.trim.fastq_trim(FASTQ(fileobj=fq), linker_3='ACGT', min_len=10, quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
@foo
aaccggttccttggaa
+
;;;;;;;;;;;;;;;;
@bar - note offset
tgtgatagctacgact
+
;;;;;;;;;;;;;;;;
@baz
tgtagatgatagatagaAGCT
+
;;;;;;;;;;;;;;;;;;;;;
''')

    def testTrim53(self):
        fq = StringIO.StringIO('''\
@foo
TTGCaaccggttccttggaaACGT
+
;;;;;;;;;;;;;;;;;;;;;;;;
@bar - note offset
aaTTGCtgtgatagctacgactACGTa
+
;;;;;;;;;;;;;;;;;;;;;;;;;;;
@baz
TAGCtgtagatgatagatagaAGCT
+
;;;;;;;;;;;;;;;;;;;;;;;;;
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.trim.fastq_trim(FASTQ(fileobj=fq), linker_5='TTGC', linker_3='ACGT', min_len=10, quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
@foo
aaccggttccttggaa
+
;;;;;;;;;;;;;;;;
@bar - note offset
tgtgatagctacgact
+
;;;;;;;;;;;;;;;;
@baz
TAGCtgtagatgatagatagaAGCT
+
;;;;;;;;;;;;;;;;;;;;;;;;;
''')

    def testTrimMinLen(self):
        fq = StringIO.StringIO('''\
@foo
TTGCaaccggttccttggaa
+
;;;;;;;;;;;;;;;;;;;;
@bar
TTGCtgtgatagctacgactaaaacc
+
;;;;;;;;;;;;;;;;;;;;;;;;;;
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.trim.fastq_trim(FASTQ(fileobj=fq), linker_5='TTGC', min_len=20, quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
@bar
tgtgatagctacgactaaaacc
+
;;;;;;;;;;;;;;;;;;;;;;
''')

    def testTrimMinOffset(self):
        fq = StringIO.StringIO('''\
@foo
aaaaaTTGCaaccggttccttggaa
+
;;;;;;;;;;;;;;;;;;;;;;;;;
@bar
TTGCtgtgatagctacgactaaaacc
+
;;;;;;;;;;;;;;;;;;;;;;;;;;
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.trim.fastq_trim(FASTQ(fileobj=fq), linker_5='TTGC', min_len=20, quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
@foo
aaaaaTTGCaaccggttccttggaa
+
;;;;;;;;;;;;;;;;;;;;;;;;;
@bar
tgtgatagctacgactaaaacc
+
;;;;;;;;;;;;;;;;;;;;;;
''')

    def testTrimPctIdentity(self):
        fq = StringIO.StringIO('''\
@foo
TTGCCaaccggttccttagaa
+
;;;;;;;;;;;;;;;;;;;;;
@bar
TTGGCtgtgatagctacgactaaaacc
+
;;;;;;;;;;;;;;;;;;;;;;;;;;;
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.trim.fastq_trim(FASTQ(fileobj=fq), linker_5='TTGGC', min_len=10, pct_identity=1.0, quiet=True, out=out)

        self.assertEqual(out.getvalue(), '''\
@foo
TTGCCaaccggttccttagaa
+
;;;;;;;;;;;;;;;;;;;;;
@bar
tgtgatagctacgactaaaacc
+
;;;;;;;;;;;;;;;;;;;;;;
''')

    def testTrimPctIdentity0_8(self):
        fq = StringIO.StringIO('''\
@foo
TTGCCaaccggttccttagaa
+
;;;;;;;;;;;;;;;;;;;;;
@bar
TTGGCtgtgatagctacgactaaaacc
+
;;;;;;;;;;;;;;;;;;;;;;;;;;;
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.trim.fastq_trim(FASTQ(fileobj=fq), linker_5='TTGGC', min_len=10, pct_identity=0.75, quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
@foo
aaccggttccttagaa
+
;;;;;;;;;;;;;;;;
@bar
tgtgatagctacgactaaaacc
+
;;;;;;;;;;;;;;;;;;;;;;
''')

    def testTrimCS(self):
        # Note: colorspace trimming should really only be applicable to the 3' end
        fq = StringIO.StringIO('''\
@foo
T012301231231013101112231
+
;;;;;;;;;;;;;;;;;;;;;;;;
@bar
T012301231231013101112111
+
;;;;;;;;;;;;;;;;;;;;;;;;
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.trim.fastq_trim(FASTQ(fileobj=fq), linker_3='112231', min_len=10, quiet=True, out=out)
        self.assertEqual(out.getvalue(), '''\
@foo
T012301231231013101
+
;;;;;;;;;;;;;;;;;;
@bar
T012301231231013101112111
+
;;;;;;;;;;;;;;;;;;;;;;;;
''')

if __name__ == '__main__':
    unittest.main()
