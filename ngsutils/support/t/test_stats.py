#!/usr/bin/env python
'''
Tests for ngsutils support / docutils
'''

import unittest
import doctest

import ngsutils.support.stats


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(ngsutils.support.stats))
    return tests

if __name__ == '__main__':
    unittest.main()
