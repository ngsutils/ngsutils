#!/usr/bin/env python
'''
Tests for bamutils cleancigar
'''

import unittest
from ngsutils.bam.t import MockRead

import ngsutils.bam
import ngsutils.bam.cleancigar


class CleanCigarTest(unittest.TestCase):

    def testCleanCigar(self):
        read = MockRead('foo', 'AAAATTTTCCCGGG', 'AAAABBBBBBBCCC', tid=1, pos=1, cigar='10M0I4M', tags=[('XA', 1)])

        code = ngsutils.bam.cleancigar.read_cleancigar(read)
        self.assertEqual(code, True)

        self.assertEqual(read.qname, 'foo')  # update occurs in place
        self.assertEqual(ngsutils.bam.cigar_tostr(read.cigar), '14M')

    def testNoClean(self):
        read = MockRead('foo', 'AAAATTTTCCCGGG', 'AAAABBBBBBBCCC', tid=1, pos=1, cigar='14M', tags=[('XA', 1)])

        code = ngsutils.bam.cleancigar.read_cleancigar(read)
        self.assertEqual(code, False)
        self.assertEqual(ngsutils.bam.cigar_tostr(read.cigar), '14M')

    def testUnmapped(self):
        read = MockRead('foo', 'AAAATTTTCCCGGG', 'AAAABBBBBBBCCC', tags=[('XA', 1)])

        code = ngsutils.bam.cleancigar.read_cleancigar(read)
        self.assertEqual(code, False)


if __name__ == '__main__':
    unittest.main()
