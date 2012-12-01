#!/usr/bin/env python
'''
Tests for gtfutils / add_xref
'''

import os
import unittest
import StringIO

import ngsutils.gtf.add_xref

gtf = os.path.join(os.path.dirname(__file__), 'test-iso.gtf')
xref = os.path.join(os.path.dirname(__file__), 'test-xref.txt')


class GTFAddXrefTest(unittest.TestCase):
    def testAddXref(self):
        out = StringIO.StringIO('')
        ngsutils.gtf.add_xref.gtf_add_xref(gtf, xref, out=out, quiet=True)
        self.assertEquals(out.getvalue().replace('\t', '|'), '''\
chr1|test|gene|1001|1100|0|+|.|gene_id "foo1"; transcript_id "txpt1"; gene_name "gene1";
chr1|test|gene|1051|1100|0|+|.|gene_id "foo1"; transcript_id "txpt2"; gene_name "gene1";
chr1|test|gene|2001|2100|0|+|.|gene_id "foo2"; transcript_id "txpt3"; gene_name "gene2";
chr1|test|gene|3001|3100|0|+|.|gene_id "foo3"; transcript_id "txpt4"; gene_name "gene3";
''')

if __name__ == '__main__':
    unittest.main()
