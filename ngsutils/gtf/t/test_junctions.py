#!/usr/bin/env python
'''
Tests for gtfutils / junctions
'''

import os
import unittest
import StringIO

import ngsutils.gtf.junctions
from ngsutils.gtf import GTF

# >test1
#          1         2         3         4         5         6         7         8         9         100
# 1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
# aaaaaaaaaCCCCCCCATGCtttttttttGCGCTTTGATCcccccccccCTGAGGGGGGGGGGGGGATCGgggggggggACTgggggggTCGAGGGGGGG

# exons:
# 10,20
# 30,40
# 50,70
# 90,100

# opt: 80-82

fa = os.path.join(os.path.dirname(__file__), 'test-junc.fa')


class GTFJunctionsTest(unittest.TestCase):
    def testJunctionsSimple(self):
        gtf = GTF(fileobj=StringIO.StringIO('''\
test1|test|exon|10|20|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
test1|test|exon|30|40|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
test1|test|exon|50|70|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
test1|test|exon|90|100|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
'''.replace('|', '\t')), quiet=True)

        valid = '''\
>test1:16-20,29-33
ATGCGCGC
>test1:16-20,49-53
ATGCCTGA
>test1:16-20,89-93
ATGCTCGA
>test1:36-40,49-53
GATCCTGA
>test1:36-40,89-93
GATCTCGA
>test1:66-70,89-93
ATCGTCGA
'''
        out = StringIO.StringIO('')
        ngsutils.gtf.junctions.gtf_junctions(gtf, fa, fragment_size=4, min_size=8, out=out, quiet=True)
        self.assertEqual(out.getvalue(), valid)

    def testJunctionsMultiExon(self):
        gtf = GTF(fileobj=StringIO.StringIO('''\
test1|test|exon|30|40|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
test1|test|exon|50|70|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
test1|test|exon|80|82|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
test1|test|exon|90|100|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
'''.replace('|', '\t')), quiet=True)

        valid = '''\
>test1:36-40,49-53
GATCCTGA
>test1:36-40,79-82,89-93
GATCACTTCGA
>test1:36-40,89-93
GATCTCGA
>test1:66-70,79-82,89-93
ATCGACTTCGA
>test1:66-70,89-93
ATCGTCGA
'''
        out = StringIO.StringIO('')
        ngsutils.gtf.junctions.gtf_junctions(gtf, fa, fragment_size=4, min_size=8, out=out, quiet=True)

        self.assertEqual(out.getvalue(), valid)

    def testJunctionsIsoforms(self):
        gtf = GTF(fileobj=StringIO.StringIO('''\
test1|test|exon|10|20|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
test1|test|exon|30|40|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
test1|test|exon|90|100|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
test1|test|exon|10|20|0|+|.|gene_id "foo1"; transcript_id "bar2"; isoform_id "iso1"
test1|test|exon|50|70|0|+|.|gene_id "foo1"; transcript_id "bar2"; isoform_id "iso1"
test1|test|exon|90|100|0|+|.|gene_id "foo1"; transcript_id "bar2"; isoform_id "iso1"
'''.replace('|', '\t')), quiet=True)

        valid = '''\
>test1:16-20,29-33
ATGCGCGC
>test1:16-20,49-53
ATGCCTGA
>test1:16-20,89-93
ATGCTCGA
>test1:36-40,49-53
GATCCTGA
>test1:36-40,89-93
GATCTCGA
>test1:66-70,89-93
ATCGTCGA
'''
        out = StringIO.StringIO('')
        ngsutils.gtf.junctions.gtf_junctions(gtf, fa, fragment_size=4, min_size=8, out=out, quiet=True)
        self.assertEqual(out.getvalue(), valid)

    def testJunctionsIsoformsKnown(self):
        gtf = GTF(fileobj=StringIO.StringIO('''\
test1|test|exon|10|20|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
test1|test|exon|30|40|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
test1|test|exon|90|100|0|+|.|gene_id "foo1"; transcript_id "bar1"; isoform_id "iso1"
test1|test|exon|10|20|0|+|.|gene_id "foo1"; transcript_id "bar2"; isoform_id "iso1"
test1|test|exon|50|70|0|+|.|gene_id "foo1"; transcript_id "bar2"; isoform_id "iso1"
test1|test|exon|90|100|0|+|.|gene_id "foo1"; transcript_id "bar2"; isoform_id "iso1"
'''.replace('|', '\t')), quiet=True)

        valid = '''\
>test1:16-20,29-33
ATGCGCGC
>test1:36-40,89-93
GATCTCGA
>test1:16-20,49-53
ATGCCTGA
>test1:66-70,89-93
ATCGTCGA
'''
        out = StringIO.StringIO('')
        ngsutils.gtf.junctions.gtf_junctions(gtf, fa, fragment_size=4, min_size=8, known=True, out=out, quiet=True)

        self.assertEqual(out.getvalue(), valid)


if __name__ == '__main__':
    unittest.main()
