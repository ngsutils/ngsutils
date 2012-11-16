#!/usr/bin/env python
'''
Tests for bedutils tofasta
'''

import os
import unittest

import ngsutils.bed.tofasta
import StringIO
from ngsutils.bed import BedFile

bedtest = BedFile(fileobj=StringIO.StringIO('''\
test1|0|10|foo1|10|+
test1|10|20|foo2|10|-
test1|0|5|foo1|10|+
test3|0|50|foo1|10|+
'''.replace('|', '\t')))

fasta = os.path.join(os.path.dirname(__file__), 'test.fa')


class FASTATest(unittest.TestCase):
    def testBedFASTA(self):

        valid = '''\
>test1:0-10
aaaaaaaaaa
>test1:10-20
cccccccccc
'''

        sio = StringIO.StringIO("")
        ngsutils.bed.tofasta.bed_tofasta(bedtest, fasta, min_size=10, stranded=False, out=sio)
        self.assertEqual(valid, sio.getvalue())
        sio.close()

    def testBedFASTAStranded(self):

        valid = '''\
>test1:0-10[+]
aaaaaaaaaa
>test1:10-20[-]
gggggggggg
'''

        sio = StringIO.StringIO("")
        ngsutils.bed.tofasta.bed_tofasta(bedtest, fasta, min_size=10, stranded=True, out=sio)
        self.assertEqual(valid, sio.getvalue())
        sio.close()


if __name__ == '__main__':
    unittest.main()
