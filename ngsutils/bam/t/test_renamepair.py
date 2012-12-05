#!/usr/bin/env python
'''
Tests for bamutils renamepair
'''

import unittest
from ngsutils.bam.t import MockRead, assertIn, assertNotIn

import ngsutils.bam.renamepair


class RenamePairTest(unittest.TestCase):
    def testRename(self):
        read = MockRead('foo/1')

        ngsutils.bam.renamepair.read_renamepair(read, '/')

        self.assertEqual(read.qname, 'foo')  # update occurs in place
        assertIn(('ZN', '1'), read.tags)  # added tag

    def testNoRename(self):
        read = MockRead('foo1')

        ngsutils.bam.renamepair.read_renamepair(read, '/')

        self.assertEqual(read.qname, 'foo1')
        assertNotIn(('ZN', '1'), read.tags)


if __name__ == '__main__':
    unittest.main()
