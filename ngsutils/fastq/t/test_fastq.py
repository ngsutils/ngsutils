#!/usr/bin/env python
'''
Tests for fastqutils / docutils
'''

import unittest
import doctest

import ngsutils.fastq
import ngsutils.fastq.fromfasta


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(ngsutils.fastq))
    tests.addTests(doctest.DocTestSuite(ngsutils.fastq.fromfasta))
    return tests

if __name__ == '__main__':
    unittest.main()
