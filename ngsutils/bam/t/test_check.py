#!/usr/bin/env python
'''
Tests for bamutils check
'''

import os
import unittest

import ngsutils.bam
import ngsutils.bam.check


class CheckTest(unittest.TestCase):
    def testCheck(self):
        ret = ngsutils.bam.check.bam_check(os.path.join(os.path.dirname(__file__), 'test.bam'), quiet=True)
        self.assertEqual(ret, True)

    def testCheckFail(self):
        ret = ngsutils.bam.check.bam_check(os.path.join(os.path.dirname(__file__), 'test1.bam'), quiet=True)
        self.assertEqual(ret, False)

if __name__ == '__main__':
    unittest.main()
