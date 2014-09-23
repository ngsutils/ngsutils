#!/usr/bin/env python
'''
Tests for fastqutils filter
'''

import unittest
import StringIO

import ngsutils.fastq.filter
from ngsutils.fastq import FASTQ


class FilterTest(unittest.TestCase):
    def testFilterNull(self):
        fq = StringIO.StringIO('''\
@foo comment
ACGTACGT
+
;;;;;;;;
@bar
ACGTACGT
+
;;;;;;;;
''')

        out = StringIO.StringIO('')
        chain = ngsutils.fastq.filter.FASTQReader(FASTQ(fileobj=fq), verbose=False)
        ngsutils.fastq.filter.fastq_filter(chain, out=out, quiet=True)
        self.assertEqual(out.getvalue(), fq.getvalue())

    def testFilterSize(self):
        fq = StringIO.StringIO('''\
@foo comment
ACGTACGTACGTACGT
+
;;;;;;;;;;;;;;;;
@bar
ACGTACGT
+
;;;;;;;;
''')

        out = StringIO.StringIO('')
        chain = ngsutils.fastq.filter.FASTQReader(FASTQ(fileobj=fq), verbose=False)
        chain = ngsutils.fastq.filter.SizeFilter(chain, 12, verbose=False)
        ngsutils.fastq.filter.fastq_filter(chain, out=out, quiet=True)
        self.assertEqual(out.getvalue(), '''\
@foo comment
ACGTACGTACGTACGT
+
;;;;;;;;;;;;;;;;
''')

    def testFilterPaired(self):
        fq = StringIO.StringIO('''\
@foo
ACGTACGTACGTACGT
+
;;;;;;;;;;;;;;;;
@foo
ACGTACGTACGTACGT
+
;;;;;;;;;;;;;;;;
@bar
ACGTACGT
+
;;;;;;;;
@baz
ACGTACGT
+
;;;;;;;;
@baz
ACGTACGT
+
;;;;;;;;
''')

        out = StringIO.StringIO('')
        chain = ngsutils.fastq.filter.FASTQReader(FASTQ(fileobj=fq), verbose=False)
        chain = ngsutils.fastq.filter.PairedFilter(chain, verbose=False)
        ngsutils.fastq.filter.fastq_filter(chain, out=out, quiet=True)
        self.assertEqual(out.getvalue(), '''\
@foo
ACGTACGTACGTACGT
+
;;;;;;;;;;;;;;;;
@foo
ACGTACGTACGTACGT
+
;;;;;;;;;;;;;;;;
@baz
ACGTACGT
+
;;;;;;;;
@baz
ACGTACGT
+
;;;;;;;;
''')

    def testFilterQual(self):
        fq = StringIO.StringIO('''\
@foo comment
ACGTACGTACGTACGT
+
+++++++++&&&&+++
@bar
ACGTACGTA
+
+++++&&++
''')

        out = StringIO.StringIO('')
        chain = ngsutils.fastq.filter.FASTQReader(FASTQ(fileobj=fq), verbose=False)
        chain = ngsutils.fastq.filter.QualFilter(chain, 7, 5, verbose=False)
        ngsutils.fastq.filter.fastq_filter(chain, out=out, quiet=True)
        self.assertEqual(out.getvalue(), '''\
@foo comment #qual
ACGTACGTACGT
+
+++++++++&&&
@bar
ACGTACGTA
+
+++++&&++
''')

    def testFilterQualIllumina(self):
        fq = StringIO.StringIO('''\
@foo
ACGTACGTACGTACGT
+
JJJJJJJJJEEEEJJJ
@bar
ACGTACGTA
+
JJJJJEEJJ
''')

        out = StringIO.StringIO('')
        chain = ngsutils.fastq.filter.FASTQReader(FASTQ(fileobj=fq), verbose=False)
        chain = ngsutils.fastq.filter.QualFilter(chain, 7, 5, illumina=True, verbose=False)
        ngsutils.fastq.filter.fastq_filter(chain, out=out, quiet=True)
        self.assertEqual(out.getvalue(), '''\
@foo #qual
ACGTACGTACGT
+
JJJJJJJJJEEE
@bar
ACGTACGTA
+
JJJJJEEJJ
''')

    def testFilterDiscard(self):
        fq = StringIO.StringIO('''\
@foo
ACGTACGTACGTACGT
+
;;;;;;;;;;;;;;;;
@bar
ACGTACGT
+
;;;;;;;;
''')

        discarded = []

        def _discard(name):
            discarded.append(name)

        out = StringIO.StringIO('')
        chain = ngsutils.fastq.filter.FASTQReader(FASTQ(fileobj=fq), verbose=False)
        chain = ngsutils.fastq.filter.SizeFilter(chain, 12, verbose=False, discard=_discard)
        ngsutils.fastq.filter.fastq_filter(chain, out=out, quiet=True)
        self.assertEqual(discarded, ['bar'])

    def testFilterWildcard(self):
        fq = StringIO.StringIO('''\
@foo
ACGTACNNACGTATTT
+
;;;;;;;;;;;;;;;;
@bar
ACGTACGT
+
;;;;;;;;
@baz
ACGTACGNACTG
+
;;;;;;;;;;;;
@quux
ACG..CGAACTG
+
;;;;;;;;;;;;
''')

        out = StringIO.StringIO('')
        chain = ngsutils.fastq.filter.FASTQReader(FASTQ(fileobj=fq), verbose=False)
        chain = ngsutils.fastq.filter.WildcardFilter(chain, 1, verbose=False)
        ngsutils.fastq.filter.fastq_filter(chain, out=out, quiet=True)

        self.assertEqual(out.getvalue(), '''\
@bar
ACGTACGT
+
;;;;;;;;
@baz
ACGTACGNACTG
+
;;;;;;;;;;;;
''')

    def testFilterSuffixQual(self):
        fq = StringIO.StringIO('''\
@foo
ACGTACGTACGTATTT
+
;;;;;;;;;;;;!!!!
@bar
ACGTACGTA
+
;;;;;;;;!
''')

        out = StringIO.StringIO('')
        chain = ngsutils.fastq.filter.FASTQReader(FASTQ(fileobj=fq), verbose=False)
        chain = ngsutils.fastq.filter.SuffixQualFilter(chain, '!', verbose=False)
        ngsutils.fastq.filter.fastq_filter(chain, out=out, quiet=True)

        self.assertEqual(out.getvalue(), '''\
@foo #suff
ACGTACGTACGT
+
;;;;;;;;;;;;
@bar #suff
ACGTACGT
+
;;;;;;;;
''')

    def testFilterTrim(self):
        fq = StringIO.StringIO('''\
@foo
ACGTACGTACGTATTT
+
;;;;;;;;;;;;;;;;
@bar
ACGTatttACGT
+
;;;;;;;;;;;;
@baz
ACGTACGTATTC
+
;;;;;;;;;;;;
@quux
ATTTAtcgtagt
+
;;;;;;;;;;;;
''')

        out = StringIO.StringIO('')
        chain = ngsutils.fastq.filter.FASTQReader(FASTQ(fileobj=fq), verbose=False)
        chain = ngsutils.fastq.filter.TrimFilter(chain, 'ATTT', 1.0, 3, verbose=False)
        ngsutils.fastq.filter.fastq_filter(chain, out=out, quiet=True)

        self.assertEqual(out.getvalue(), '''\
@foo #trim
ACGTACGTACGT
+
;;;;;;;;;;;;
@bar #trim
ACGT
+
;;;;
@baz
ACGTACGTATTC
+
;;;;;;;;;;;;
''')

    def testFilterTrim2(self):
        '''Test to make sure that we don't trim out the middle of a read'''
        fq = StringIO.StringIO('''\
@foo
ACGTACGTACGTATTATT
+
;;;;;;;;;;;;;;;;;;
@bar
CCGGATTAGGCCGGCC
+
;;;;;;;;;;;;;;;;
@baz
CCGGATTATTATTATT
+
;;;;;;;;;;;;;;;;
''')

        out = StringIO.StringIO('')
        chain = ngsutils.fastq.filter.FASTQReader(FASTQ(fileobj=fq), verbose=False)
        chain = ngsutils.fastq.filter.TrimFilter(chain, 'ATTATT', 1.0, 4, verbose=False)
        ngsutils.fastq.filter.fastq_filter(chain, out=out, quiet=True)

        self.assertEqual(out.getvalue(), '''\
@foo #trim
ACGTACGTACGT
+
;;;;;;;;;;;;
@bar
CCGGATTAGGCCGGCC
+
;;;;;;;;;;;;;;;;
@baz #trim
CCGGATTATT
+
;;;;;;;;;;
''')

    def testFilterTrimRepeat(self):
        '''properly trim repeats'''
        fq = StringIO.StringIO('''\
@foo
ACGTACGTAAAAAAAAAA
+
;;;;;;;;;;;;;;;;;;
@bar
CCGGATTAGGCCCAAA
+
;;;;;;;;;;;;;;;;
''')

        out = StringIO.StringIO('')
        chain = ngsutils.fastq.filter.FASTQReader(FASTQ(fileobj=fq), verbose=False)
        chain = ngsutils.fastq.filter.TrimFilter(chain, 'AAAAAAAAAAAAAAAAAA', 1.0, 4, verbose=False)
        ngsutils.fastq.filter.fastq_filter(chain, out=out, quiet=True)

        self.assertEqual(out.getvalue(), '''\
@foo #trim
ACGTACGT
+
;;;;;;;;
@bar
CCGGATTAGGCCCAAA
+
;;;;;;;;;;;;;;;;
''')

    def testFilterTrimCS(self):
        fq = StringIO.StringIO('''\
@foo
T0123012301231122
+
;;;;;;;;;;;;;;;;
@bar
T012301231122
+
;;;;;;;;;;;;
@baz
T012301231102
+
;;;;;;;;;;;;
@quux
T1122Atcgtagt
+
;;;;;;;;;;;;
''')

        out = StringIO.StringIO('')
        chain = ngsutils.fastq.filter.FASTQReader(FASTQ(fileobj=fq), verbose=False)
        chain = ngsutils.fastq.filter.TrimFilter(chain, '1122', 0.8, 3, verbose=False)
        ngsutils.fastq.filter.fastq_filter(chain, out=out, quiet=True)

        self.assertEqual(out.getvalue(), '''\
@foo #trim
T012301230123
+
;;;;;;;;;;;;
@bar #trim
T01230123
+
;;;;;;;;
@baz
T012301231102
+
;;;;;;;;;;;;
''')

if __name__ == '__main__':
    unittest.main()
