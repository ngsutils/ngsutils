#!/usr/bin/env python
'''
Tests for bedutils reduce
'''

import unittest
import StringIO

import ngsutils.bed.reduce
from ngsutils.bed import BedFile
from ngsutils.bam.t import _matches

bedtest = BedFile(fileobj=StringIO.StringIO('''\
chr1|100|150|foo1|10|+
chr1|140|200|foo2|10|+
chr1|500|550|foo3|10|+
chr1|560|600|foo4|10|-
chr1|600|700|foo5|10|+
'''.replace('|', '\t')))


class ReduceTest(unittest.TestCase):
    def testReduce(self):
        valids = '''\
chr1|100|200|foo1,foo2|20|+
chr1|500|550|foo3|10|+
chr1|560|700|foo4,foo5|20|+
'''.replace('|', '\t').split('\n')

        out = StringIO.StringIO('')
        ngsutils.bed.reduce.bed_reduce(bedtest, extend=(0, 0), stranded=False, count=False, out=out)

        self.assertTrue(_matches(valids, out.getvalue().split('\n')))

    def testReduceCount(self):
        valids = '''\
chr1|100|200|foo1,foo2|2|+
chr1|500|550|foo3|1|+
chr1|560|700|foo4,foo5|2|+
'''.replace('|', '\t').split('\n')

        out = StringIO.StringIO('')
        ngsutils.bed.reduce.bed_reduce(bedtest, extend=(0, 0), stranded=False, count=True, out=out)

        self.assertTrue(_matches(valids, out.getvalue().split('\n')))

    def testReduceStranded(self):
        valids = '''\
chr1|100|200|foo1,foo2|20|+
chr1|500|550|foo3|10|+
chr1|560|600|foo4|10|-
chr1|600|700|foo5|10|+
'''.replace('|', '\t').split('\n')

        out = StringIO.StringIO('')
        ngsutils.bed.reduce.bed_reduce(bedtest, extend=(0, 0), stranded=True, count=False, out=out)
        self.assertTrue(_matches(valids, out.getvalue().split('\n')))

    def testReduceExtend(self):
        valids = '''\
chr1|90|210|foo1,foo2|20|+
chr1|490|710|foo3,foo4,foo5|30|+
'''.replace('|', '\t').split('\n')

        out = StringIO.StringIO('')
        ngsutils.bed.reduce.bed_reduce(bedtest, extend=(10, 10), stranded=False, out=out)
        self.assertTrue(_matches(valids, out.getvalue().split('\n')))

    def testReduceExtendAsync(self):
        valids = '''\
chr1|90|220|foo1,foo2|20|+
chr1|490|720|foo3,foo4,foo5|30|+
'''.replace('|', '\t').split('\n')

        out = StringIO.StringIO('')
        ngsutils.bed.reduce.bed_reduce(bedtest, extend=(10, 20), stranded=False, out=out)
        self.assertTrue(_matches(valids, out.getvalue().split('\n')))

    def testReduceExtendClip(self):
        valids = '''\
chr1|100|200|foo1,foo2|20|+
chr1|500|700|foo3,foo4,foo5|30|+
'''.replace('|', '\t').split('\n')

        out = StringIO.StringIO('')
        ngsutils.bed.reduce.bed_reduce(bedtest, extend=(10, 10), stranded=False, clip=True, out=out)
        self.assertTrue(_matches(valids, out.getvalue().split('\n')))


if __name__ == '__main__':
    unittest.main()
