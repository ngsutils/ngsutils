#!/usr/bin/env python
'''
Tests for gtfutils / docutils
'''

import unittest
import doctest
import StringIO

import ngsutils.gtf
from ngsutils.gtf import GTF


class GTFTest(unittest.TestCase):
    def testGTF(self):
        src = StringIO.StringIO('''\
chr1|test|exon|1001|1100|0|+|.|gene_id "foo"; transcript_id "bar";
chr1|test|exon|1201|1300|0|+|.|gene_id "foo"; transcript_id "bar";
chr1|test|CDS|1051|1247|0|+|1|gene_id "foo"; transcript_id "bar";
chr1|test|start_codon|1051|1053|0|+|.|gene_id "foo"; transcript_id "bar";
chr1|test|stop_codon|1248|1250|0|+|.|gene_id "foo"; transcript_id "bar";
'''.replace('|', '\t'))

        gtf = GTF(fileobj=src, quiet=True)
        for gene in gtf.genes:
            self.assertEqual(gene.gene_id, 'foo')
            for transcript in gene.transcripts:
                self.assertEqual(transcript.transcript_id, 'bar')
                self.assertEqual(list(transcript.exons), [(1000, 1100), (1200, 1300)])

    def testGeneOnly(self):
        src = StringIO.StringIO('''\
chr1|test|gene|1001|1100|0|+|.|gene_id "foo1"; transcript_id "bar1";
chr1|test|gene|2001|2100|0|+|.|gene_id "foo2"; transcript_id "bar2";
'''.replace('|', '\t'))

        gtf = GTF(fileobj=src, quiet=True)
        genes = list(gtf.genes)
        self.assertEqual(['foo1', 'foo2'], [x.gene_id for x in genes])
        # gene start / end
        self.assertEqual((1000, 1100), (genes[0].start, genes[0].end))
        t = list(genes[0].transcripts)[0]
        # apply start / end to transcript and as an exon
        self.assertEqual((1000, 1100), (t.start, t.end))
        self.assertEqual((1000, 1100), (t.exons[0]))

        t = list(genes[1].transcripts)[0]
        # apply start / end to transcript and as an exon
        self.assertEqual((2000, 2100), (t.start, t.end))
        self.assertEqual((2000, 2100), (t.exons[0]))

    def testFind(self):
        src = StringIO.StringIO('''\
chr1|test|gene|1001|1100|0|+|.|gene_id "foo1"; transcript_id "bar1";
chr1|test|gene|2001|2100|0|+|.|gene_id "foo2"; transcript_id "bar2";
chr1|test|gene|3001|3100|0|+|.|gene_id "foo3"; transcript_id "bar3";
'''.replace('|', '\t'))

        gtf = GTF(fileobj=src, quiet=True)
        genes = list(gtf.genes)
        self.assertEqual(['foo1', 'foo2', 'foo3'], [x.gene_id for x in genes])

        found = False
        for gene in gtf.find('chr1', 1900, 2200):
            found = True
            self.assertEqual(gene.gid, 'foo2')
        self.assertTrue(found)

        found = False
        for gene in gtf.find('chr1', 1900, 2000):
            found = True
            self.assertEqual(gene.gid, 'foo2')
        self.assertTrue(found)

        found = False
        for gene in gtf.find('chr1', 2050, 2200):
            found = True
            self.assertEqual(gene.gid, 'foo2')
        self.assertTrue(found)

    def testGTFNoIso(self):
        src = StringIO.StringIO('''\
chr1|test|exon|1001|1100|0|+|.|gene_id "foo"; transcript_id "bar1";
chr1|test|exon|1201|1300|0|+|.|gene_id "foo"; transcript_id "bar1";
chr1|test|exon|1401|1500|0|+|.|gene_id "foo"; transcript_id "bar1";
chr1|test|CDS|1051|1447|0|+|1|gene_id "foo"; transcript_id "bar1";
chr1|test|start_codon|1051|1053|0|+|.|gene_id "foo"; transcript_id "bar1";
chr1|test|stop_codon|1448|1450|0|+|.|gene_id "foo"; transcript_id "bar1";
chr1|test|exon|1001|1100|0|+|.|gene_id "foo"; transcript_id "bar2";
chr1|test|exon|1401|1500|0|+|.|gene_id "foo"; transcript_id "bar2";
chr1|test|CDS|1051|1447|0|+|1|gene_id "foo"; transcript_id "bar2";
chr1|test|start_codon|1051|1053|0|+|.|gene_id "foo"; transcript_id "bar2";
chr1|test|stop_codon|1448|1450|0|+|.|gene_id "foo"; transcript_id "bar2";
'''.replace('|', '\t'))

        gtf = GTF(fileobj=src, quiet=True)
        genes = list(gtf.genes)
        self.assertEqual('foo', genes[0].gene_id)
        transcripts = list(genes[0].transcripts)
        self.assertEqual(len(transcripts), 2)
        self.assertEqual(list(genes[0].regions), [(1, 1000, 1100, True, 'bar1,bar2'), (2, 1200, 1300, False, 'bar1'), (3, 1400, 1500, True, 'bar1,bar2')])

    def testGTFIso(self):
        src = StringIO.StringIO('''\
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
'''.replace('|', '\t'))

        gtf = GTF(fileobj=src, quiet=True)
        genes = list(gtf.genes)
        self.assertEqual(len(genes), 1)
        self.assertEqual(['foo1'], [g.gene_id for g in genes])
        self.assertEqual(['iso1'], [g.gid for g in genes])  # gid is the actual unique id (should be constant for isoforms)
        transcripts = list(genes[0].transcripts)
        self.assertEqual(['bar1', 'bar2'], [t.transcript_id for t in transcripts])
        self.assertEqual(len(transcripts), 2)
        self.assertEqual(list(genes[0].regions), [(1, 1000, 1100, True, 'bar1,bar2'), (2, 1200, 1300, False, 'bar1'), (3, 1400, 1500, True, 'bar1,bar2')])


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite(ngsutils.gtf))
    return tests

if __name__ == '__main__':
    unittest.main()
