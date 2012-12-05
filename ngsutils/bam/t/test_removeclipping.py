#!/usr/bin/env python
'''
Tests for bamutils removeclipping
'''

import unittest
from ngsutils.bam.t import MockRead, assertIn

import ngsutils.bam
import ngsutils.bam.removeclipping


class RemoveClippingTest(unittest.TestCase):
    def testUnmapped(self):
        read = MockRead('foo', 'AAAATTTTCCCGGG', 'AAAABBBBBBBCCC', tags=[('XA', 1)])

        code = ngsutils.bam.removeclipping.read_removeclipping(read)
        self.assertEqual(code, 1)
        self.assertEqual(read.qname, 'foo')  # update occurs in place
        assertIn(('XA', 1), read.tags)  # added tag
        self.assertEqual(read.seq, 'AAAATTTTCCCGGG')
        self.assertEqual(read.qual, 'AAAABBBBBBBCCC')

    def testRemoveClipping(self):
        read = MockRead('foo', 'AAAATTTTCCCGGG', 'AAAABBBBBBBCCC', tags=[('XA', 1)], tid=1, pos=1, cigar='4S7M3S')

        code = ngsutils.bam.removeclipping.read_removeclipping(read)
        self.assertEqual(code, 2)
        self.assertEqual(read.qname, 'foo')  # update occurs in place
        assertIn(('XA', 1), read.tags)  # added tag
        assertIn(('ZA', 4), read.tags)  # added tag
        assertIn(('ZB', 3), read.tags)  # added tag
        assertIn(('ZC', 0.5), read.tags)  # added tag
        self.assertEqual(read.seq, 'TTTTCCC')
        self.assertEqual(read.qual, 'BBBBBBB')

    def testNoClipping(self):
        read = MockRead('foo', 'AAAATTTTCCCGGG', 'AAAABBBBBBBCCC', tags=[('XA', 1)], tid=1, pos=1, cigar='14M')

        code = ngsutils.bam.removeclipping.read_removeclipping(read)
        self.assertEqual(code, 0)
        self.assertEqual(read.qname, 'foo')  # update occurs in place
        assertIn(('XA', 1), read.tags)  # added tag
        self.assertEqual(read.seq, 'AAAATTTTCCCGGG')
        self.assertEqual(read.qual, 'AAAABBBBBBBCCC')


if __name__ == '__main__':
    unittest.main()
