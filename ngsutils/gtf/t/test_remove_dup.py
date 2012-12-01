#!/usr/bin/env python
'''
Tests for gtfutils / removedup
'''

import os
import unittest
import StringIO

import ngsutils.gtf.remove_dup

fname = os.path.join(os.path.dirname(__file__), 'test1.gtf')


class GTFRemoveDupTest(unittest.TestCase):
    def testRemoveDup(self):
        out = StringIO.StringIO('')
        retval = ngsutils.gtf.remove_dup.gtf_remove_dup(fname, out=out, quiet=True)
        self.assertEquals((6, 6), retval)
        self.assertEquals(out.getvalue().replace('\t', '|'), '''\
chr1|test|exon|1001|1100|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
chr1|test|exon|1201|1300|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
chr1|test|exon|1401|1500|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
chr1|test|CDS|1051|1447|0|+|1|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
chr1|test|start_codon|1051|1053|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
chr1|test|stop_codon|1448|1450|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
''')

if __name__ == '__main__':
    unittest.main()
