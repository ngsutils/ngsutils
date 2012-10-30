#!/usr/bin/env python
'''
Tests for bamutils tag
'''

import unittest
from mock import *

import ngsutils.bam
import ngsutils.bam.removeclipping


class RemoveClippingTest(unittest.TestCase):
    def testUnmapped(self):
        read = Mock()
        read.qname = 'foo'
        read.tags = [('XA', 1)]
        read.seq = 'AAAATTTTCCCGGG'
        read.qual = 'AAAABBBBBBBCCC'
        read.is_unmapped = True

        code = ngsutils.bam.removeclipping.read_removeclipping(read)
        self.assertEqual(code, 1)

        self.assertEqual(read.qname, 'foo')  # update occurs in place
        self.assertIn(('XA', 1), read.tags)  # added tag
        self.assertEqual(read.seq, 'AAAATTTTCCCGGG')
        self.assertEqual(read.qual, 'AAAABBBBBBBCCC')

    def testRemoveClipping(self):
        read = Mock()
        read.qname = 'foo'
        read.tags = [('XA', 1)]
        read.seq = 'AAAATTTTCCCGGG'
        read.qual = 'AAAABBBBBBBCCC'
        read.is_unmapped = False
        read.cigar = ngsutils.bam.cigar_fromstr('4S7M3S')

        code = ngsutils.bam.removeclipping.read_removeclipping(read)
        self.assertEqual(code, 2)

        self.assertEqual(read.qname, 'foo')  # update occurs in place
        self.assertIn(('XA', 1), read.tags)  # added tag
        self.assertIn(('ZA', 4), read.tags)  # added tag
        self.assertIn(('ZB', 3), read.tags)  # added tag
        self.assertIn(('ZC', 0.5), read.tags)  # added tag
        self.assertEqual(read.seq, 'TTTTCCC')
        self.assertEqual(read.qual, 'BBBBBBB')

    def testNoClipping(self):
        read = Mock()
        read.qname = 'foo'
        read.tags = [('XA', 1)]
        read.seq = 'AAAATTTTCCCGGG'
        read.qual = 'AAAABBBBBBBCCC'
        read.is_unmapped = False
        read.cigar = ngsutils.bam.cigar_fromstr('14M')

        code = ngsutils.bam.removeclipping.read_removeclipping(read)
        self.assertEqual(code, 0)

        self.assertEqual(read.qname, 'foo')  # update occurs in place
        self.assertIn(('XA', 1), read.tags)  # added tag
        self.assertEqual(read.seq, 'AAAATTTTCCCGGG')
        self.assertEqual(read.qual, 'AAAABBBBBBBCCC')


if __name__ == '__main__':
    unittest.main()
