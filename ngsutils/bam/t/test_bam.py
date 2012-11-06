#!/usr/bin/env python
'''
Tests for bamutils / docutils
'''

import unittest
import doctest

import ngsutils.bam
import ngsutils.bam.convertregion


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(ngsutils.bam))
    tests.addTests(doctest.DocTestSuite(ngsutils.bam.convertregion))
    return tests

if __name__ == '__main__':
    unittest.main()
