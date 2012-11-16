#!/usr/bin/env python
'''
Tests for bedutils refcount
'''

import unittest
import StringIO

import ngsutils.bed.refcount
from ngsutils.bed import BedFile

bedtest1 = BedFile('test1', fileobj=StringIO.StringIO('''\
chr1|100|150|foo1|10|+
chr1|140|200|foo2|10|+
chr1|500|550|foo3|10|+
chr1|560|600|foo4|10|-
chr1|600|700|foo5|10|+
'''.replace('|', '\t')))

bedtest2 = BedFile('test2', fileobj=StringIO.StringIO('''\
chr1|90|120|bar1|10|+
chr1|140|200|bar2|10|+
chr1|510|520|bar3|10|+
chr1|500|660|bar4|10|-
chr1|520|570|bar5|10|+
'''.replace('|', '\t')))

bedtest3 = BedFile('test3', fileobj=StringIO.StringIO('''\
chr1|100|150|baz1|10|+
chr1|750|850|baz2|10|+
'''.replace('|', '\t')))

bedref = BedFile(fileobj=StringIO.StringIO('''\
chr1|100|200|foo1|10|+
chr1|500|700|foo2|10|-
'''.replace('|', '\t')))


class RefCountTest(unittest.TestCase):
    def testRefCount(self):
        valid = '''\
#chrom|start|end|name|score|strand|test1|test2|test3|present
chr1|100|200|foo1|10|+|2|2|1|3
chr1|500|700|foo2|10|-|3|3|0|2
'''.replace('|', '\t')

        out = StringIO.StringIO('')
        ngsutils.bed.refcount.bed_refcount(bedref, [bedtest1, bedtest2, bedtest3], stranded=False, out=out)

        self.assertEqual(valid, out.getvalue())

    def testRefCountStranded(self):
        valid = '''\
#chrom|start|end|name|score|strand|test1|test2|test3|present
chr1|100|200|foo1|10|+|2|2|1|3
chr1|500|700|foo2|10|-|1|1|0|2
'''.replace('|', '\t')

        out = StringIO.StringIO('')
        ngsutils.bed.refcount.bed_refcount(bedref, [bedtest1, bedtest2, bedtest3], stranded=True, out=out)

        self.assertEqual(valid, out.getvalue())


if __name__ == '__main__':
    unittest.main()
