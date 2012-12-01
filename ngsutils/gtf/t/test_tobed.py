#!/usr/bin/env python
'''
Tests for gtfutils / tobed
'''

import unittest
import StringIO

import ngsutils.gtf.tobed
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
'''.replace('|', '\t')), quiet=True)


class GTFTest(unittest.TestCase):
    def testGeneToBed(self):
        out = StringIO.StringIO('')
        ngsutils.gtf.tobed.gtf_genes_tobed(gtf, out=out)
        self.assertEquals(out.getvalue().replace('\t', '|'), 'chr1|1000|1500|foo1|0|+\n')

    def testExonToBed(self):
        out = StringIO.StringIO('')
        ngsutils.gtf.tobed.gtf_exons_tobed(gtf, out=out)
        self.assertEquals(out.getvalue().replace('\t', '|'), '''\
chr1|1000|1100|foo1.1|0|+
chr1|1200|1300|foo1.2|0|+
chr1|1400|1500|foo1.3|0|+
''')

    def testRegionToBed(self):
        out = StringIO.StringIO('')
        ngsutils.gtf.tobed.gtf_regions_tobed(gtf, out=out)
        self.assertEquals(out.getvalue().replace('\t', '|'), '''\
chr1|1000|1100|foo1.const.1|2|+
chr1|1200|1300|foo1.alt.2|1|+
chr1|1400|1500|foo1.const.3|2|+
''')


if __name__ == '__main__':
    unittest.main()
