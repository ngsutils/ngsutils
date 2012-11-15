#!/usr/bin/env python
'''
Tests for bamutils / docutils
'''

import unittest
import doctest

import ngsutils.bam
import ngsutils.bam.convertregion
import ngsutils.bam.count.count


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(ngsutils.bam))
    tests.addTests(doctest.DocTestSuite(ngsutils.bam.convertregion))
    tests.addTests(doctest.DocTestSuite(ngsutils.bam.count.count))
    return tests

if __name__ == '__main__':
    unittest.main()
