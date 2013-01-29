#!/usr/bin/env python
'''
Tests for fastqutils split
'''

import os
import unittest

import ngsutils.fastq.split
from ngsutils.fastq import FASTQ


class SplitTest(unittest.TestCase):
    def testSplit(self):
        fname = os.path.join(os.path.dirname(__file__), 'test.fastq')
        templ = os.path.join(os.path.dirname(__file__), 'test_templ')

        ngsutils.fastq.split.fastq_split(fname, templ, 2, quiet=True)

        self.assertTrue(os.path.exists('%s.1.fastq' % templ))
        self.assertTrue(os.path.exists('%s.2.fastq' % templ))

        fq1 = FASTQ('%s.1.fastq' % templ)
        fq2 = FASTQ('%s.2.fastq' % templ)

        names1 = [x.fullname for x in fq1.fetch(quiet=True)]
        self.assertEqual(names1, ['foo /1', 'foo /2', 'baz /1', 'baz /2'])

        names2 = [x.fullname for x in fq2.fetch(quiet=True)]
        self.assertEqual(names2, ['bar /1', 'bar /2'])

        fq1.close()
        fq2.close()
        os.unlink('%s.1.fastq' % templ)
        os.unlink('%s.2.fastq' % templ)

    def testSplitGz(self):
        fname = os.path.join(os.path.dirname(__file__), 'test.fastq')
        templ = os.path.join(os.path.dirname(__file__), 'test_templ')

        ngsutils.fastq.split.fastq_split(fname, templ, 2, gz=True, quiet=True)

        self.assertTrue(os.path.exists('%s.1.fastq.gz' % templ))
        self.assertTrue(os.path.exists('%s.2.fastq.gz' % templ))

        os.unlink('%s.1.fastq.gz' % templ)
        os.unlink('%s.2.fastq.gz' % templ)

    def testSplitUnpaired(self):
        fname = os.path.join(os.path.dirname(__file__), 'test.fastq')
        templ = os.path.join(os.path.dirname(__file__), 'test_templ')

        ngsutils.fastq.split.fastq_split(fname, templ, 2, ignore_pairs=True, quiet=True)

        self.assertTrue(os.path.exists('%s.1.fastq' % templ))
        self.assertTrue(os.path.exists('%s.2.fastq' % templ))

        fq1 = FASTQ('%s.1.fastq' % templ)
        fq2 = FASTQ('%s.2.fastq' % templ)

        names1 = [x.name for x in fq1.fetch(quiet=True)]
        self.assertEqual(names1, ['foo', 'bar', 'baz'])

        names2 = [x.name for x in fq2.fetch(quiet=True)]
        self.assertEqual(names2, ['foo', 'bar', 'baz'])

        fq1.close()
        fq2.close()
        os.unlink('%s.1.fastq' % templ)
        os.unlink('%s.2.fastq' % templ)

    def testSplitThree(self):
        fname = os.path.join(os.path.dirname(__file__), 'test.fastq')
        templ = os.path.join(os.path.dirname(__file__), 'test_templ')

        ngsutils.fastq.split.fastq_split(fname, templ, 3, ignore_pairs=True, quiet=True)

        self.assertTrue(os.path.exists('%s.1.fastq' % templ))
        self.assertTrue(os.path.exists('%s.2.fastq' % templ))
        self.assertTrue(os.path.exists('%s.3.fastq' % templ))

        fq1 = FASTQ('%s.1.fastq' % templ)
        fq2 = FASTQ('%s.2.fastq' % templ)
        fq3 = FASTQ('%s.3.fastq' % templ)

        names1 = [x.fullname for x in fq1.fetch(quiet=True)]
        self.assertEqual(names1, ['foo /1', 'bar /2'])

        names2 = [x.fullname for x in fq2.fetch(quiet=True)]
        self.assertEqual(names2, ['foo /2', 'baz /1'])

        names3 = [x.fullname for x in fq3.fetch(quiet=True)]
        self.assertEqual(names3, ['bar /1', 'baz /2'])

        fq1.close()
        fq2.close()
        fq3.close()

        os.unlink('%s.1.fastq' % templ)
        os.unlink('%s.2.fastq' % templ)
        os.unlink('%s.3.fastq' % templ)

if __name__ == '__main__':
    unittest.main()
