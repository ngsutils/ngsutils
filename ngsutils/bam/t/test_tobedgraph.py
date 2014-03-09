#!/usr/bin/env python
'''
Tests for bamutils tobedgraph
'''

import os
import unittest

import ngsutils.bam
import ngsutils.bam.tobedgraph
import StringIO


class ToBEDGraphTest(unittest.TestCase):
    def setUp(self):
        self.bam = ngsutils.bam.bam_open(os.path.join(os.path.dirname(__file__), 'test.bam'))

    def tearDown(self):
        self.bam.close()

    def testBEDGraph(self):
        sio = StringIO.StringIO("")

        ngsutils.bam.tobedgraph.bam_tobedgraph(self.bam, out=sio)
        self.assertEqual(sio.getvalue(), '''\
chr1\t99\t149\t1
chr1\t174\t199\t2
chr1\t399\t499\t1
chr1\t699\t724\t2
chr1\t724\t774\t1
''')
        sio.close()

    def testBEDGraphPlus(self):
        sio = StringIO.StringIO("")

        ngsutils.bam.tobedgraph.bam_tobedgraph(self.bam, strand='+', out=sio)
        self.assertEqual(sio.getvalue(), '''\
chr1\t99\t149\t1
chr1\t174\t199\t2
chr1\t399\t499\t1
chr1\t699\t724\t2
''')
        sio.close()

    def testBEDGraphMinus(self):
        sio = StringIO.StringIO("")

        ngsutils.bam.tobedgraph.bam_tobedgraph(self.bam, strand='-', out=sio)
        self.assertEqual(sio.getvalue(), '''\
chr1\t724\t774\t1
''')
        sio.close()

if __name__ == '__main__':
    unittest.main()
