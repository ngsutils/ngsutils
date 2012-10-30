#!/usr/bin/env python
'''
Tests for bamutils tag
'''

import unittest
from mock import *

import ngsutils.bam
import ngsutils.bam.cleancigar


class CleanCigarTest(unittest.TestCase):

    def testCleanCigar(self):
        read = Mock()
        read.qname = 'foo'
        read.tags = [('XA', 1)]
        read.seq = 'AAAATTTTCCCGGG'
        read.qual = 'AAAABBBBBBBCCC'
        read.is_unmapped = False
        read.cigar = ngsutils.bam.cigar_fromstr('10M0I4M')

        code = ngsutils.bam.cleancigar.read_cleancigar(read)
        self.assertEqual(code, True)

        self.assertEqual(read.qname, 'foo')  # update occurs in place
        self.assertEqual(ngsutils.bam.cigar_tostr(read.cigar), '14M')

    def testNoClean(self):
        read = Mock()
        read.qname = 'foo'
        read.tags = [('XA', 1)]
        read.seq = 'AAAATTTTCCCGGG'
        read.qual = 'AAAABBBBBBBCCC'
        read.is_unmapped = False
        read.cigar = ngsutils.bam.cigar_fromstr('14M')

        code = ngsutils.bam.cleancigar.read_cleancigar(read)
        self.assertEqual(code, False)
        self.assertEqual(ngsutils.bam.cigar_tostr(read.cigar), '14M')

    def testUnmapped(self):
        read = Mock()
        read.qname = 'foo'
        read.tags = [('XA', 1)]
        read.seq = 'AAAATTTTCCCGGG'
        read.qual = 'AAAABBBBBBBCCC'
        read.is_unmapped = True

        code = ngsutils.bam.cleancigar.read_cleancigar(read)
        self.assertEqual(code, False)


if __name__ == '__main__':
    unittest.main()
