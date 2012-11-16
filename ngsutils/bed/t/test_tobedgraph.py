#!/usr/bin/env python
'''
Tests for bedutils tobedgraph
'''

import unittest

import ngsutils.bed.tobedgraph
import StringIO
from ngsutils.bed import BedFile

bedtest = BedFile(fileobj=StringIO.StringIO('''\
test1|10|20|foo1|10|+
test1|10|20|foo2|10|-
test1|15|25|foo3|10|+
test1|100|150|foo4|1|+
'''.replace('|', '\t')))


class BedGraphTest(unittest.TestCase):
    def testBedGraph(self):
        valid = '''\
test1|10|15|2
test1|15|20|3
test1|20|25|1
test1|100|150|1
'''.replace('|', '\t')
        sio = StringIO.StringIO("")
        ngsutils.bed.tobedgraph.bed_tobedgraph(bedtest, out=sio)
        self.assertEqual(valid, sio.getvalue())
        sio.close()

    def testBedGraphNormalize(self):
        valid = '''\
test1|10|15|4
test1|15|20|6
test1|20|25|2
test1|100|150|2
'''.replace('|', '\t')
        sio = StringIO.StringIO("")
        ngsutils.bed.tobedgraph.bed_tobedgraph(bedtest, out=sio, normalize=2)
        self.assertEqual(valid, sio.getvalue())
        sio.close()

    def testBedGraphPlus(self):
        valid = '''\
test1|10|15|1
test1|15|20|2
test1|20|25|1
test1|100|150|1
'''.replace('|', '\t')
        sio = StringIO.StringIO("")
        ngsutils.bed.tobedgraph.bed_tobedgraph(bedtest, out=sio, only_strand='+')
        self.assertEqual(valid, sio.getvalue())
        sio.close()

    def testBedGraphMinus(self):
        valid = '''\
test1|10|20|1
'''.replace('|', '\t')
        sio = StringIO.StringIO("")
        ngsutils.bed.tobedgraph.bed_tobedgraph(bedtest, out=sio, only_strand='-')
        self.assertEqual(valid, sio.getvalue())
        sio.close()

if __name__ == '__main__':
    unittest.main()
