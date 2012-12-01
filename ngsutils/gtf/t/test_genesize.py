#!/usr/bin/env python
'''
Tests for gtfutils / genesize
'''

import unittest
import StringIO

import ngsutils.gtf.genesize
from ngsutils.gtf import GTF

gtf = GTF(fileobj=StringIO.StringIO('''\
chr1|test|exon|1001|1100|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
chr1|test|exon|1201|1300|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
chr1|test|exon|1401|1500|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
chr1|test|CDS|1051|1447|0|+|1|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
chr1|test|start_codon|1051|1053|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
chr1|test|stop_codon|1448|1450|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
chr1|test|exon|1001|1100|0|+|.|gene_id "foo2"; transcript_id "bar2"; isoform_id "iso1"
chr1|test|exon|1401|1500|0|+|.|gene_id "foo2"; transcript_id "bar2"; isoform_id "iso1"
chr1|test|CDS|1051|1447|0|+|1|gene_id "foo2"; transcript_id "bar2"; isoform_id "iso1"
chr1|test|start_codon|1051|1053|0|+|.|gene_id "foo2"; transcript_id "bar2"; isoform_id "iso1"
chr1|test|stop_codon|1448|1450|0|+|.|gene_id "foo2"; transcript_id "bar2"; isoform_id "iso1"
chr1|test|exon|2001|2100|0|+|.|gene_id "foo3"; transcript_id "bar1"; isoform_id "iso2"
chr1|test|exon|2201|2300|0|+|.|gene_id "foo3"; transcript_id "bar1"; isoform_id "iso2"
chr1|test|exon|2401|2500|0|+|.|gene_id "foo3"; transcript_id "bar1"; isoform_id "iso2"
chr1|test|CDS|2051|2447|0|+|1|gene_id "foo3"; transcript_id "bar1"; isoform_id "iso2"
chr1|test|start_codon|2051|2053|0|+|.|gene_id "foo3"; transcript_id "bar1"; isoform_id "iso2"
chr1|test|stop_codon|2448|2450|0|+|.|gene_id "foo3"; transcript_id "bar1"; isoform_id "iso2"
'''.replace('|', '\t')), quiet=True)


class GTFGeneSizeTest(unittest.TestCase):
    def testGeneSize(self):
        out = StringIO.StringIO('')
        ngsutils.gtf.genesize.gtf_genesize(gtf, out=out)
        self.assertEquals(out.getvalue().replace('\t', '|'), '''\
#gene|genomic-size|transcript-size
foo1|500|300
foo3|500|300
''')


if __name__ == '__main__':
    unittest.main()
