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

		chain = _TestChain([read,])
		chain = ngsutils.bam.tag.Suffix(chain, 'bar')  # add 'bar' as a suffix to all reads

		for r in chain.filter():
			self.assertEqual(r.qname, 'foobar')

	def testTag(self):
		read = Mock()
		read.qname = 'foo'
		read.tags=[('ZZ', 'val')]

		chain = _TestChain([read,])
		chain = ngsutils.bam.tag.Tag(chain, 'ZY:Z:value1')
		chain = ngsutils.bam.tag.Tag(chain, 'ZX:value2')

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
