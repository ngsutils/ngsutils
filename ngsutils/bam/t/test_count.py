#!/usr/bin/env python
'''
Tests for bamutils count
'''

import unittest
import StringIO

import ngsutils.bam
import ngsutils.bam.count
import ngsutils.bam.count.models

from ngsutils.bam.t import MockBam

testbam1 = MockBam(['chr1'])
testbam1.add_read('foo1', 'A' * 50, tid=0, pos=100, cigar='50M', tags=[('IH', 2)])
testbam1.add_read('foo1', 'A' * 50, tid=0, pos=1000, cigar='50M', tags=[('IH', 2)])
testbam1.add_read('foo2', 'A' * 50, tid=0, pos=101, cigar='50M')
testbam1.add_read('foo3', 'A' * 50, tid=0, pos=101, cigar='50M')
testbam1.add_read('foo4', 'A' * 50, tid=0, pos=103, cigar='50M')
testbam1.add_read('foo5', 'A' * 50, tid=0, pos=200, cigar='25M')
testbam1.add_read('foo6', 'A' * 50, tid=0, pos=200, cigar='25M')
testbam1.add_read('foo7', 'A' * 50, tid=0, pos=225, cigar='25M')
testbam1.add_read('foo8', 'A' * 50, tid=0, pos=225, cigar='25M')
testbam1.add_read('foo9', 'A' * 50, tid=0, pos=300, cigar='50M')
testbam1.add_read('foo10', 'A' * 50, tid=0, pos=301, cigar='50M')
testbam1.add_read('foo11', 'A' * 50, tid=0, pos=302, cigar='50M', is_reverse=True)
testbam1.add_read('foo12', 'A' * 50, tid=0, pos=303, cigar='50M', is_reverse=True)


class CountTest(unittest.TestCase):
    def testCountBed(self):
        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|100|150|foo|1|+
chr1|300|350|foo|1|-
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, out=out, quiet=True)

        valid = '''\
## input
## model bed *fileobj*
## stranded False
## multiple complete
#chrom|start|end|name|score|strand|length|count
chr1|100|150|foo|1|+|50|4
chr1|300|350|foo|1|-|50|4
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)

    def testCountBedCoverage(self):
        # first, look only in 20 base range where the first 4 reads map
        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|110|130|foo|1|+
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, coverage=True, out=out, quiet=True)

        valid = '''\
## input
## model bed *fileobj*
## stranded False
## multiple complete
#chrom|start|end|name|score|strand|length|count|coverage mean|coverage stdev|coverage median
chr1|110|130|foo|1|+|20|4|4.0|0.0|4
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)

    def testCountBedCoverage2(self):
        # look between 199-299, where there should be 2 reads coverage (4 reads)

        # note: there may still be an issue with an off-by one error for the region...
        # I'm not sure if it is in the bedfile iterator, or testing basecaller/pileup used

        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|200|249|foo|1|+
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, coverage=True, out=out, quiet=True)

        valid = '''\
## input
## model bed *fileobj*
## stranded False
## multiple complete
#chrom|start|end|name|score|strand|length|count|coverage mean|coverage stdev|coverage median
chr1|200|249|foo|1|+|49|4|2.0|0.0|2.0
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)

    def testCountBedNormAll(self):
        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|100|150|foo|1|+
chr1|300|350|foo|1|-
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, norm='all', out=out, quiet=True)

        valid = '''\
## input
## model bed *fileobj*
## stranded False
## multiple complete
## norm all 12.0
## CPM-factor 1.2e-05
#chrom|start|end|name|score|strand|length|count|count (CPM)
chr1|100|150|foo|1|+|50|4|333333.333333
chr1|300|350|foo|1|-|50|4|333333.333333
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)

    def testCountBedNormMedian(self):
        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|100|150|foo|1|+
chr1|300|350|foo|1|-
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, norm='median', out=out, quiet=True)

        valid = '''\
## input
## model bed *fileobj*
## stranded False
## multiple complete
## norm median 4.0
## CPM-factor 4e-06
#chrom|start|end|name|score|strand|length|count|count (CPM)
chr1|100|150|foo|1|+|50|4|1000000.0
chr1|300|350|foo|1|-|50|4|1000000.0
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)

    def testCountBedNormMapped(self):
        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|100|150|foo|1|+
chr1|300|350|foo|1|-
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, norm='mapped', out=out, quiet=True)

        valid = '''\
## input
## model bed *fileobj*
## stranded False
## multiple complete
## norm mapped 8.0
## CPM-factor 8e-06
#chrom|start|end|name|score|strand|length|count|count (CPM)
chr1|100|150|foo|1|+|50|4|500000.0
chr1|300|350|foo|1|-|50|4|500000.0
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)

    def testCountBedNormQuartile(self):
        #Not tested - mock bam is too small - tests for this in doctest
        pass

    def testCountBedUniq(self):
        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|100|150|foo|1|+
chr1|300|350|foo|1|-
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, uniq_only=True, out=out, quiet=True)

        valid = '''\
## input
## model bed *fileobj*
## stranded False
## multiple complete
#chrom|start|end|name|score|strand|length|count
chr1|100|150|foo|1|+|50|3
chr1|300|350|foo|1|-|50|4
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)
    def testCountBedUniq(self):
        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|100|150|foo|1|+
chr1|300|350|foo|1|-
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, uniq_only=True, out=out, quiet=True)

        valid = '''\
## input
## model bed *fileobj*
## stranded False
## multiple complete
#chrom|start|end|name|score|strand|length|count
chr1|100|150|foo|1|+|50|3
chr1|300|350|foo|1|-|50|4
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)

    def testCountBedStranded(self):
        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|100|150|foo|1|+
chr1|300|350|foo|1|-
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, stranded=True, out=out, quiet=True)

        valid = '''\
## input
## model bed *fileobj*
## stranded True
## multiple complete
#chrom|start|end|name|score|strand|length|count
chr1|100|150|foo|1|+|50|4
chr1|300|350|foo|1|-|50|2
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)

    def testCountBedMultiIgnore(self):
        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|100|150|foo|1|+
chr1|300|350|foo|1|-
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, multiple='ignore', out=out, quiet=True)

        valid = '''\
## input
## model bed *fileobj*
## stranded False
## multiple ignore
#chrom|start|end|name|score|strand|length|count
chr1|100|150|foo|1|+|50|3
chr1|300|350|foo|1|-|50|4
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)

    def testCountBedMultiPartial(self):
        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|100|150|foo|1|+
chr1|300|350|foo|1|-
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, multiple='partial', out=out, quiet=True)

        valid = '''\
## input
## model bed *fileobj*
## stranded False
## multiple partial
#chrom|start|end|name|score|strand|length|count
chr1|100|150|foo|1|+|50|3.5
chr1|300|350|foo|1|-|50|4
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)

    def testCountBedRPKM(self):
        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|100|150|foo|1|+
chr1|300|350|foo|1|-
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, rpkm=True, norm="mapped", out=out, quiet=True)

        valid = '''\
## input
## model bed *fileobj*
## stranded False
## multiple complete
## norm mapped 8.0
## CPM-factor 8e-06
#chrom|start|end|name|score|strand|length|count|count (CPM)|RPKM
chr1|100|150|foo|1|+|50|4|500000.0|10000000.0
chr1|300|350|foo|1|-|50|4|500000.0|10000000.0
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)

    def testCountBedWhitelist(self):
        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|100|150|foo|1|+
chr1|300|350|foo|1|-
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, out=out, quiet=True, whitelist=['foo2', 'foo3', 'foo9', 'foo10', 'foo11'])

        valid = '''\
## input
## model bed *fileobj*
## stranded False
## multiple complete
#chrom|start|end|name|score|strand|length|count
chr1|100|150|foo|1|+|50|2
chr1|300|350|foo|1|-|50|3
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)

    def testCountBedBlacklist(self):
        counter = ngsutils.bam.count.models['bed'](fileobj=StringIO.StringIO('''
chr1|100|150|foo|1|+
chr1|300|350|foo|1|-
'''.replace('|', '\t')))

        out = StringIO.StringIO('')
        counter.count(testbam1, out=out, quiet=True, blacklist=['foo2', 'foo3', 'foo9', 'foo10', 'foo11'])

        valid = '''\
## input
## model bed *fileobj*
## stranded False
## multiple complete
#chrom|start|end|name|score|strand|length|count
chr1|100|150|foo|1|+|50|2
chr1|300|350|foo|1|-|50|1
'''.replace('|', '\t')

        self.assertEquals(out.getvalue(), valid)


def dump(s, t):
    print 'valid:'
    print s.replace('\t', '|')
    print 'output:'
    print t.replace('\t', '|')

if __name__ == '__main__':
    unittest.main()
