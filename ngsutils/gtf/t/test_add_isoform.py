#!/usr/bin/env python
'''
Tests for gtfutils / add_isoform
'''

import os
import unittest
import StringIO

import ngsutils.gtf.add_isoform

gtf = os.path.join(os.path.dirname(__file__), 'test-iso.gtf')
iso = os.path.join(os.path.dirname(__file__), 'test-iso.txt')


class GTFAddIsoformTest(unittest.TestCase):
    def testAddIsoform(self):
        out = StringIO.StringIO('')
        ngsutils.gtf.add_isoform.gtf_add_isoform(gtf, iso, out=out, quiet=True)
        self.assertEquals(out.getvalue().replace('\t', '|'), '''\
chr1|test|gene|1001|1100|0|+|.|gene_id "foo1"; transcript_id "txpt1"; isoform_id "1";
chr1|test|gene|1051|1100|0|+|.|gene_id "foo1"; transcript_id "txpt2"; isoform_id "1";
chr1|test|gene|2001|2100|0|+|.|gene_id "foo2"; transcript_id "txpt3"; isoform_id "2";
chr1|test|gene|3001|3100|0|+|.|gene_id "foo3"; transcript_id "txpt4"; isoform_id "3";
''')

if __name__ == '__main__':
    unittest.main()
