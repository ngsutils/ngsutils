#!/usr/bin/env python
'''
Tests for bamutils tag
'''

import unittest
from mock import *

import ngsutils.bam.tag


class TagTest(unittest.TestCase):
    def testSuffix(self):
        read = Mock()
        read.qname = 'foo'

        chain = _TestChain([read, ])
        chain = ngsutils.bam.tag.Suffix(chain, 'bar')  # add 'bar' as a suffix to all reads

        for r in chain.filter():
            self.assertEqual(r.qname, 'foobar')

    def testCufflinksXS(self):
        read1 = Mock()
        read1.qname = 'foo'
        read1.is_unmapped = False
        read1.is_reverse = True
        read1.tags = [('ZZ', 1)]

        read2 = Mock()
        read2.qname = 'bar'
        read2.is_unmapped = False
        read2.is_reverse = False
        read2.tags = [('ZZ', 1)]

        read3 = Mock()
        read3.qname = 'baz'
        read3.is_unmapped = True
        read3.tags = [('ZZ', 1)]

        chain = _TestChain([read1, read2, read3, ])
        chain = ngsutils.bam.tag.CufflinksXS(chain)

        for r in chain.filter():
            self.assertIn(('ZZ', 1), r.tags)
            if r.qname == 'foo':
                self.assertIn(('XS', '-'), r.tags)
            elif r.qname == 'bar':
                self.assertIn(('XS', '+'), r.tags)
            elif r.qname == 'baz':
                self.assertNotIn(('XS', '+'), r.tags)
                self.assertNotIn(('XS', '-'), r.tags)

    def testTag(self):
        read = Mock()
        read.qname = 'foo'
        read.tags = [('ZZ', 'val')]

        chain = _TestChain([read, ])
        chain = ngsutils.bam.tag.Tag(chain, 'ZY:Z:value1')  # add ZY:Z (should remove the :Z)
        chain = ngsutils.bam.tag.Tag(chain, 'ZX:value2')  # add ZX

        for r in chain.filter():
            found_orig = False
            found_1 = False
            found_2 = False

            for k, v in r.tags:
                if k == 'ZZ' and v == 'val':
                    found_orig = True
                if k == 'ZY' and v == 'value1':
                    found_1 = True
                if k == 'ZX' and v == 'value2':
                    found_2 = True

            self.assertTrue(found_orig)
            self.assertTrue(found_1)
            self.assertTrue(found_2)


class _TestChain(object):
    def __init__(self, vals):
        self.vals = vals

    def filter(self):
        for v in self.vals:
            yield v

if __name__ == '__main__':
    unittest.main()
