#!/usr/bin/env python
'''
Tests for bamutils tag
'''

import unittest
from mock import *

import ngsutils.bam.renamepair


class RenamePairTest(unittest.TestCase):
    def testRename(self):
        read = Mock()
        read.qname = 'foo/1'
        read.tags = []

        ngsutils.bam.renamepair.read_renamepair(read, '/')

        self.assertEqual(read.qname, 'foo')  # update occurs in place
        self.assertIn(('ZN', '1'), read.tags)  # added tag


if __name__ == '__main__':
    unittest.main()
