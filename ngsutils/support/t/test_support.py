#!/usr/bin/env python
'''
Tests for ngsutils support / docutils
'''

import unittest
import doctest

import ngsutils.support.ngs_utils


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(ngsutils.support))
    tests.addTests(doctest.DocTestSuite(ngsutils.support.ngs_utils))
    return tests


class CountsTest(unittest.TestCase):
    def testBins(self):
        counts = ngsutils.support.Counts()
        counts.add(1)
        counts.add(2)
        counts.add(3)
        counts.add(4)
        counts.add(5)
        self.assertTrue(counts.max(), 5)
        self.assertTrue(counts.mean(), 3)

    def testBins2(self):
        counts = ngsutils.support.Counts()
        counts.add(1)
        counts.add(1)
        counts.add(1)
        counts.add(1)
        counts.add(1)
        counts.add(1)
        counts.add(2)
        counts.add(3)
        counts.add(4)
        counts.add(5)

        self.assertTrue(counts.mean(), 2)

if __name__ == '__main__':
    unittest.main()
