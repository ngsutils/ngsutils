#!/usr/bin/env python
'''
Tests for bedutils frombasecall
'''

import unittest
import StringIO

import ngsutils.bed.frombasecall


class FromBasecallTest(unittest.TestCase):
    def testBasecall(self):
        valid = 'chr1|10|11\nchr1|11|12\n'.replace('|', '\t')
        basecall = StringIO.StringIO('''\
chr1|11|A|50|A|T|1|0.5|1|25|0|0|25|0|0|0
chr1|12|C|50|C|G|1|0.5|1|0|25|25|0|0|0|0
'''.replace('|', '\t'))

        out = StringIO.StringIO('')
        ngsutils.bed.frombasecall._bed_frombasecall(basecall, out=out)
        self.assertEqual(out.getvalue(), valid)


if __name__ == '__main__':
    unittest.main()
