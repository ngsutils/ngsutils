#!/usr/bin/env python
'''
Tests for barcode_split

Note: These tests use small barcodes, so some of the test sequences are pretty
arbitrary. In real-life, the barcodes should be much longer (12-16bp).

'''

import unittest
import os

import ngsutils.fastq.barcode_split
from ngsutils.support import FASTA
from ngsutils.fastq import FASTQ


barcodes = {
    'tag1': ('ATAT', '5'),
    'tag2': ('TGTG', '5'),
    'tag3': ('CTCT', '3')
}

barcodes2 = {
    'tag1': ('AATTAA', '5'),
    'tag2': ('GGTTCC', '5'),
    'tag3': ('CCAACC', '3')
}


class BarcodeSplitTest(unittest.TestCase):
    def test_check_tags_5(self):
        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'ATATaaaatttt', 0, 0, False)
        self.assertTrue(valid)
        self.assertEqual(results[0], 'tag1')
        self.assertTrue(results[2])

        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'TGTGaaaatttt', 0, 0, False)
        self.assertTrue(valid)
        self.assertEqual(results[0], 'tag2')
        self.assertTrue(results[2])

    def test_check_tags_5_multi(self):
        "pulls out only the 5' ATAT, so the next ATAT is kept in-tact"
        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'ATATATATaaaatttt', 0, 0, False)
        self.assertTrue(valid)
        self.assertEqual(results[0], 'tag1')
        self.assertTrue(results[2])

    def test_check_tags_3(self):
        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'aattttaaCTCT', 0, 0, False)
        # print valid, results
        self.assertTrue(valid)
        self.assertEqual(results[0], 'tag3')
        self.assertTrue(results[2])

    def test_check_tags_5_fail(self):
        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'aaaattttATAT', 0, 0, False)
        self.assertFalse(valid)
        self.assertEqual(results[0], 'tag1')  # not at right location
        self.assertTrue(results[2])

    def test_check_tags_5_revcomp(self):
        # matches ATAT in rev-comp
        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'aaaattttATAT', 0, 0, True)

        self.assertTrue(valid)
        self.assertEqual(results[0], 'tag1')
        self.assertFalse(results[2])

    def test_check_tags_5_mm(self):
        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'ATTTacgtacgt', 0, 0, False)
        self.assertFalse(valid)
        self.assertEqual(results[0], 'tag1')
        self.assertTrue(results[2])

        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'ATTTacgtacgt', 1, 0, False)
        self.assertTrue(valid)
        self.assertEqual(results[0], 'tag1')
        self.assertTrue(results[2])

    def test_check_tags_5_pos(self):
        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'gATATacgtacgt', 0, 0, False)
        self.assertFalse(valid)
        self.assertEqual(results[0], 'tag1')
        self.assertTrue(results[2])

        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'gATATacgtacgt', 0, 1, False)
        self.assertTrue(valid)
        self.assertEqual(results[0], 'tag1')
        self.assertTrue(results[2])

    def test_check_tags_3_mm(self):
        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'ttttaaaaggggCGCT', 0, 0, False)
        self.assertFalse(valid)
        self.assertEqual(results[0], 'tag3')
        self.assertTrue(results[2])

        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'ttttaaaaggggCGCT', 1, 0, False)
        self.assertTrue(valid)
        self.assertEqual(results[0], 'tag3')
        self.assertTrue(results[2])

    def test_check_tags_3_pos(self):
        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'ttttaaaaggggCTCTg', 0, 0, False)
        self.assertFalse(valid)
        self.assertEqual(results[0], 'tag3')
        self.assertTrue(results[2])

        valid, results = ngsutils.fastq.barcode_split.check_tags(barcodes, 'ttttaaaaggggCTCTg', 0, 1, False)
        self.assertTrue(valid)
        self.assertEqual(results[0], 'tag3')
        self.assertTrue(results[2])

    def test_splitFasta(self):
        path = os.path.dirname(__file__)
        ngsutils.fastq.barcode_split.fastx_barcode_split(FASTA(os.path.join(path, 'test_barcodes.fasta')), os.path.join(path, 'out.%s.fasta'), barcodes2)
        self.assert_fasta_contains(os.path.join(path, 'out.%s.fasta'), {
            'missing': 'foo-rc foo1 bar1 baz1 foo2 bar2 baz2 foo1-rc foo2-rc',
            'tag1': 'foo',
            'tag2': 'bar',
            'tag3': 'baz'
            })
        self._unlink_fastx(os.path.join(path, 'out.%s.fasta'), 'missing tag1 tag2 tag3'.split())

    def test_splitFastaRevComp(self):
        path = os.path.dirname(__file__)
        ngsutils.fastq.barcode_split.fastx_barcode_split(FASTA(os.path.join(path, 'test_barcodes.fasta')), os.path.join(path, 'out.%s.fasta'), barcodes2, allow_revcomp=True)
        self.assert_fasta_contains(os.path.join(path, 'out.%s.fasta'), {
            'missing': 'foo1 bar1 baz1 foo2 bar2 baz2 foo1-rc foo2-rc',
            'tag1': 'foo foo-rc',
            'tag2': 'bar',
            'tag3': 'baz'
            })
        self._unlink_fastx(os.path.join(path, 'out.%s.fasta'), 'missing tag1 tag2 tag3'.split())

    def test_splitFastaEdit(self):
        path = os.path.dirname(__file__)
        ngsutils.fastq.barcode_split.fastx_barcode_split(FASTA(os.path.join(path, 'test_barcodes.fasta')), os.path.join(path, 'out.%s.fasta'), barcodes2, allow_revcomp=True, edits=1)
        self.assert_fasta_contains(os.path.join(path, 'out.%s.fasta'), {
            'missing': 'foo2 bar2 baz2 foo2-rc',
            'tag1': 'foo foo-rc foo1 foo1-rc',
            'tag2': 'bar bar1',
            'tag3': 'baz baz1'
            })
        self._unlink_fastx(os.path.join(path, 'out.%s.fasta'), 'missing tag1 tag2 tag3'.split())

    def test_splitFastaOffset(self):
        path = os.path.dirname(__file__)
        ngsutils.fastq.barcode_split.fastx_barcode_split(FASTA(os.path.join(path, 'test_barcodes.fasta')), os.path.join(path, 'out.%s.fasta'), barcodes2, allow_revcomp=True, pos=1)
        self.assert_fasta_contains(os.path.join(path, 'out.%s.fasta'), {
            'missing': 'foo1 bar1 baz1 foo1-rc',
            'tag1': 'foo foo-rc foo2 foo2-rc',
            'tag2': 'bar bar2',
            'tag3': 'baz baz2'
            })
        self._unlink_fastx(os.path.join(path, 'out.%s.fasta'), 'missing tag1 tag2 tag3'.split())

    def test_splitFastq(self):
        path = os.path.dirname(__file__)
        ngsutils.fastq.barcode_split.fastx_barcode_split(FASTQ(os.path.join(path, 'test_barcodes.fastq')), os.path.join(path, 'out.%s.fastq'), barcodes2, allow_revcomp=True)
        self.assert_fastq_contains(os.path.join(path, 'out.%s.fastq'), {
            'missing': ('quux', '', ''),
            'tag1': ('foo foo-rc', 'atcgatcgatcgatcg atcgatcgatcgatcg', 'AAAAAAAAAAAAAAAA AAAAAAAAAAAAAAAA'),
            'tag2': ('bar', 'gctagctagctagcta', 'AAAAAAAAAAAAAAAA'),
            'tag3': ('baz', 'acgtacgtacgtacgt', 'AAAAAAAAAAAAAAAA')
            })
        self._unlink_fastx(os.path.join(path, 'out.%s.fastq'), 'missing tag1 tag2 tag3'.split())

    def _unlink_fastx(self, base, names):
        for name in names:
            os.unlink(base % name)

    def assert_fasta_contains(self, base, args):
        for tag in args:
            valid = args[tag].split()
            fa = FASTA(base % tag)
            count = 0
            for read in fa.fetch():
                if read.name in valid:
                    count += 1
                else:
                    self.assertEqual('extra read in %s' % tag, read.name)

            self.assertEqual(count, len(valid))

    def assert_fastq_contains(self, base, args):
        for tag in args:
            valid = args[tag][0].split()
            seq_qual = {}
            if args[tag][1]:
                for n, s, q in zip(valid, args[tag][1].split(), args[tag][2].split()):
                    seq_qual[n] = (s, q)

            fq = FASTQ(base % tag)
            count = 0
            for read in fq.fetch():
                if read.name in valid:
                    count += 1
                    if seq_qual:
                        self.assertEqual(seq_qual[read.name], (read.seq, read.qual))
                else:
                    self.assertEqual('extra read in %s' % tag, read.name)

            self.assertEqual(count, len(valid))



if __name__ == '__main__':
    unittest.main()
