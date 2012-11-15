#!/usr/bin/env python
'''
Tests for bedutils fromprimers

This test does not run automatically because it relies on a connection
to the UCSC servers to run the query. It would be possible to mock this
out, but at this point it is probably overkill.
'''

import unittest
import StringIO

import ngsutils.bed.fromprimers


class FromPrimersTest(unittest.TestCase):
    def testPrimersTab(self):
        valids = '''\
chr22|34304505|34304952|1|0|+
chr22|34304505|34304952|2|0|-
'''.replace('|', '\t').split('\n')
        tabinput = StringIO.StringIO('''\
TACAGATTGATGATGCATGAAATGGG|CATGAGTGGCTCCTAAAGCAGCTGC
CATGAGTGGCTCCTAAAGCAGCTGC|TACAGATTGATGATGCATGAAATGGG
'''.replace('|', '\t').upper())

        out = StringIO.StringIO('')
        ngsutils.bed.fromprimers._insilico_pcr_tab(tabinput, 'hg19', out=out, quiet=True)
        self.assertEqual(valids, out.getvalue().split('\n'))

    def testPrimersTabFlip(self):
        valid = 'chr22|34304505|34304952|1|0|+\n'.replace('|', '\t')
        tabinput = StringIO.StringIO('''\
TACAGATTGATGATGCATGAAATGGG|GCAGCTGCTTTAGGAGCCACTCATG
'''.replace('|', '\t').upper())

        out = StringIO.StringIO('')
        ngsutils.bed.fromprimers._insilico_pcr_tab(tabinput, 'hg19', flip=True, out=out, quiet=True)
        self.assertEqual(out.getvalue(), valid)

    def testPrimersTabMissing(self):
        valid = ''.replace('|', '\t')
        tabinput = StringIO.StringIO('''\
atcgatcgatcgatcgatcg|tgcatgcatgcatgcatgac
'''.replace('|', '\t').upper())

        out = StringIO.StringIO('')
        ngsutils.bed.fromprimers._insilico_pcr_tab(tabinput, 'hg19', out=out, quiet=True)
        self.assertEqual(out.getvalue(), valid)

    def testPrimersFasta(self):
        valids = 'chr22|34304505|34304952|foo/bar|0|+\n'.replace('|', '\t')
        tabinput = StringIO.StringIO('''\
>foo/bar/1
TACAGATTGATGATGCATGAAATGGG
>foo/bar/2
CATGAGTGGCTCCTAAAGCAGCTGC
''')

        out = StringIO.StringIO('')
        ngsutils.bed.fromprimers._insilico_pcr_fasta(tabinput, 'hg19', out=out, quiet=True)
        self.assertEqual(valids, out.getvalue())


if __name__ == '__main__':
    unittest.main()
