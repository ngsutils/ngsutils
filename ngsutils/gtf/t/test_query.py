#!/usr/bin/env python
'''
Tests for gtfutils / query

- Slightly redundant with fetch test in test_gtf.py, but that's okay :)
'''

import os
import unittest

import ngsutils.gtf.query
from ngsutils.gtf import GTF

fname = os.path.join(os.path.dirname(__file__), 'test1.gtf')


class GTFQueryTest(unittest.TestCase):
    def testQuery(self):
        genes = list(ngsutils.gtf.query.gtf_query(GTF(fname, cache_enabled=False), 'chr1', 1000, 2000))
        self.assertEquals(str(genes[0]), 'foo1(iso1) chr1:1000-2500[+]')

if __name__ == '__main__':
    unittest.main()
