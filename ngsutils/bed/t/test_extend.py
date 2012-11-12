#!/usr/bin/env python
'''
Tests for bedutils extend
'''

import unittest
import StringIO

from ngsutils.bed import BedFile
import ngsutils.bed.extend


testbed = BedFile(fileobj=StringIO.StringIO('''
chr1|10|90|foo|1|+
chr1|10|90|foo|1|-
chr1|100|150|foo|1|+
chr1|200|250|foo|1|-
'''.replace('|', '\t')))


class ExtendTest(unittest.TestCase):
    def testAbsolute(self):
        valid = '''chr1|10|110|foo|1|+
chr1|0|90|foo|1|-
chr1|100|200|foo|1|+
chr1|150|250|foo|1|-
'''.replace('|', '\t')

        out = StringIO.StringIO('')
        ngsutils.bed.extend.bed_extend(testbed, 100, relative=False, out=out)
        self.assertEqual(out.getvalue(), valid)

    def testRelative(self):
        valid = '''chr1|10|190|foo|1|+
chr1|0|90|foo|1|-
chr1|100|250|foo|1|+
chr1|100|250|foo|1|-
'''.replace('|', '\t')

        out = StringIO.StringIO('')
        ngsutils.bed.extend.bed_extend(testbed, 100, relative=True, out=out)
        self.assertEqual(out.getvalue(), valid)

if __name__ == '__main__':
    unittest.main()
