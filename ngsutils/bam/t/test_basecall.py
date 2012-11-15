#!/usr/bin/env python
'''
Tests for bamutils basecall
'''

import os
import StringIO
import unittest

from ngsutils.bam.t import MockBam
from ngsutils.bed import BedFile
import ngsutils.bam.basecall


class BaseCallTest(unittest.TestCase):
    def testBaseCall(self):
        bam = MockBam(['test2'])
        bam.add_read('foo1', 'atcgatcg', '........', 0, 0, cigar='8M')
        bam.add_read('foo2', 'atcgatcg', 'AAAAAAAA', 0, 4, cigar='8M')
        bam.add_read('foo3', 'accgatcg', '########', 0, 4, cigar='8M', tags=[('IH', 2), ('HI', 1)])

        out = StringIO.StringIO('')
        ngsutils.bam.basecall.bam_basecall(bam, os.path.join(os.path.dirname(__file__), 'test.fa'), out=out)

        valid = '''chrom|pos|ref|count|consensus call|minor call|ave mappings|entropy|A|C|G|T|N|Deletions|Gaps|Insertions|Inserts
test2|1|A|1|A||1.0|0.91686730441|1|0|0|0|0|0|0|0|
test2|2|T|1|T||1.0|3.25488750216|0|0|0|1|0|0|0|0|
test2|3|C|1|C||1.0|3.25488750216|0|1|0|0|0|0|0|0|
test2|4|G|1|G||1.0|0.91686730441|0|0|1|0|0|0|0|0|
test2|5|A|3|A||1.33333333333|1.65545182665|3|0|0|0|0|0|0|0|
test2|6|T|3|T|C|1.33333333333|2.9097148334|0|1|0|2|0|0|0|0|
test2|7|C|3|C||1.33333333333|4.85048518563|0|3|0|0|0|0|0|0|
test2|8|G|3|G||1.33333333333|1.65545182665|0|0|3|0|0|0|0|0|
test2|9|A|2|A||1.5|1.36179536943|2|0|0|0|0|0|0|0|
test2|10|T|2|T||1.5|4.24102241591|0|0|0|2|0|0|0|0|
test2|11|C|2|C||1.5|4.24102241591|0|2|0|0|0|0|0|0|
test2|12|G|2|G||1.5|1.36179536943|0|0|2|0|0|0|0|0|
'''.replace('|', '\t')
        self.assertEqual(valid, out.getvalue())

    def testBaseCallInDel(self):
        bam = MockBam(['test2'])
        bam.add_read('foo1', 'atcgaTtcg', '.........', 0, 0, cigar='5M1I3M')  # inserted T @6
        bam.add_read('foo2', 'atcgacg', 'AAAAAAA', 0, 4, cigar='5M1D2M')  # deleted T @ 10

        out = StringIO.StringIO('')
        ngsutils.bam.basecall.bam_basecall(bam, os.path.join(os.path.dirname(__file__), 'test.fa'), out=out)

        valid = '''chrom|pos|ref|count|consensus call|minor call|ave mappings|entropy|A|C|G|T|N|Deletions|Gaps|Insertions|Inserts
test2|1|A|1|A||1.0|0.91686730441|1|0|0|0|0|0|0|0|
test2|2|T|1|T||1.0|3.25488750216|0|0|0|1|0|0|0|0|
test2|3|C|1|C||1.0|3.25488750216|0|1|0|0|0|0|0|0|
test2|4|G|1|G||1.0|0.91686730441|0|0|1|0|0|0|0|0|
test2|5|A|2|A||1.0|1.36179536943|2|0|0|0|0|0|0|0|
test2|6|T|2|T||1.0|4.24102241591|0|0|0|2|0|0|0|1|T:1
test2|7|C|2|C||1.0|4.24102241591|0|2|0|0|0|0|0|0|
test2|8|G|2|G||1.0|1.36179536943|0|0|2|0|0|0|0|0|
test2|9|A|1|A||1.0|0.91686730441|1|0|0|0|0|0|0|0|
test2|10|T|0|N||1.0|0|0|0|0|0|0|1|0|0|
test2|11|C|1|C||1.0|3.25488750216|0|1|0|0|0|0|0|0|
test2|12|G|1|G||1.0|0.91686730441|0|0|1|0|0|0|0|0|
'''.replace('|', '\t')
        self.assertEqual(valid, out.getvalue())

    def testBaseCallMinQual(self):
        bam = MockBam(['test2'])
        bam.add_read('foo1', 'atcgatcg', '........', 0, 0, cigar='8M')
        bam.add_read('foo2', 'atcgatcg', 'AAAAAAAA', 0, 4, cigar='8M')
        bam.add_read('foo3', 'accgatcg', '########', 0, 4, cigar='8M', tags=[('IH', 2), ('HI', 1)])

        out = StringIO.StringIO('')
        ngsutils.bam.basecall.bam_basecall(bam, os.path.join(os.path.dirname(__file__), 'test.fa'), min_qual=10, out=out)

        valid = '''chrom|pos|ref|count|consensus call|minor call|ave mappings|entropy|A|C|G|T|N|Deletions|Gaps|Insertions|Inserts
test2|1|A|1|A||1.0|0.91686730441|1|0|0|0|0|0|0|0|
test2|2|T|1|T||1.0|3.25488750216|0|0|0|1|0|0|0|0|
test2|3|C|1|C||1.0|3.25488750216|0|1|0|0|0|0|0|0|
test2|4|G|1|G||1.0|0.91686730441|0|0|1|0|0|0|0|0|
test2|5|A|2|A||1.0|1.36179536943|2|0|0|0|0|0|0|0|
test2|6|T|2|T||1.0|4.24102241591|0|0|0|2|0|0|0|0|
test2|7|C|2|C||1.0|4.24102241591|0|2|0|0|0|0|0|0|
test2|8|G|2|G||1.0|1.36179536943|0|0|2|0|0|0|0|0|
test2|9|A|1|A||1.0|0.91686730441|1|0|0|0|0|0|0|0|
test2|10|T|1|T||1.0|3.25488750216|0|0|0|1|0|0|0|0|
test2|11|C|1|C||1.0|3.25488750216|0|1|0|0|0|0|0|0|
test2|12|G|1|G||1.0|0.91686730441|0|0|1|0|0|0|0|0|
'''.replace('|', '\t')

        self.assertEqual(valid, out.getvalue())

    def testBaseCallStrand(self):
        bam = MockBam(['test2'])
        bam.add_read('foo1', 'atcgatcg', '........', 0, 0, cigar='8M')
        bam.add_read('foo2', 'atcgatcg', 'AAAAAAAA', 0, 4, cigar='8M')
        bam.add_read('foo3', 'accgatcg', '########', 0, 4, cigar='8M', is_reverse=True)
        bam.add_read('foo4', 'atcgactgatcg', '############', 0, 0, cigar='12M')

        out = StringIO.StringIO('')
        ngsutils.bam.basecall.bam_basecall(bam, os.path.join(os.path.dirname(__file__), 'test.fa'), showstrand=True, out=out)

        valid = '''chrom|pos|ref|count|consensus call|minor call|ave mappings|entropy|A|C|G|T|N|Deletions|Gaps|Insertions|Inserts|+ strand %|A minor %|C minor %|G minor %|T minor %|N minor %|Deletion minor %|Insertion minor %
test2|1|A|2|A||1.0|1.36179536943|2|0|0|0|0|0|0|0||1.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0
test2|2|T|2|T||1.0|4.24102241591|0|0|0|2|0|0|0|0||1.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0
test2|3|C|2|C||1.0|4.24102241591|0|2|0|0|0|0|0|0||1.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0
test2|4|G|2|G||1.0|1.36179536943|0|0|2|0|0|0|0|0||1.0|0.0|0.0|0.0|0.0|0.0|0.0|0.0
test2|5|A|4|A||1.0|1.87433193885|4|0|0|0|0|0|0|0||0.75|0.25|0.0|0.0|0.0|0.0|0.0|0.0
test2|6|T|4|T/C||1.0|2.94335833285|0|2|0|2|0|0|0|0||0.75|0.0|0.5|0.0|0.0|0.0|0.0|0.0
test2|7|C|4|C|T|1.0|3.45990045591|0|3|0|1|0|0|0|0||0.75|0.0|0.333333333333|0.0|0.0|0.0|0.0|0.0
test2|8|G|4|G||1.0|1.87433193885|0|0|4|0|0|0|0|0||0.75|0.0|0.0|0.25|0.0|0.0|0.0|0.0
test2|9|A|3|A||1.0|1.65545182665|3|0|0|0|0|0|0|0||0.666666666667|0.333333333333|0.0|0.0|0.0|0.0|0.0|0.0
test2|10|T|3|T||1.0|4.85048518563|0|0|0|3|0|0|0|0||0.666666666667|0.0|0.0|0.0|0.333333333333|0.0|0.0|0.0
test2|11|C|3|C||1.0|4.85048518563|0|3|0|0|0|0|0|0||0.666666666667|0.0|0.333333333333|0.0|0.0|0.0|0.0|0.0
test2|12|G|3|G||1.0|1.65545182665|0|0|3|0|0|0|0|0||0.666666666667|0.0|0.0|0.333333333333|0.0|0.0|0.0|0.0
'''.replace('|', '\t')

        self.assertEqual(valid, out.getvalue())

    def testBaseCallRegion(self):
        bam = MockBam(['test2'])
        bam.add_read('foo1', 'atcgatcg', '........', 0, 0, cigar='8M')
        bam.add_read('foo2', 'atcgatcg', 'AAAAAAAA', 0, 4, cigar='8M')
        bam.add_read('foo3', 'accgatcg', '########', 0, 4, cigar='8M', is_reverse=True)
        bam.add_read('foo4', 'atcgactgatcg', '############', 0, 0, cigar='12M')

        bed = BedFile(fileobj=StringIO.StringIO('''\
test2|4|7
'''.replace('|', '\t')))

        #  remember: bed is 0-based, basecall is 1-based, so 4->7, corresponds to 5->8

        out = StringIO.StringIO('')
        ngsutils.bam.basecall.bam_basecall(bam, os.path.join(os.path.dirname(__file__), 'test.fa'), showstrand=True, regions=bed, out=out)

        valid = '''chrom|pos|ref|count|consensus call|minor call|ave mappings|entropy|A|C|G|T|N|Deletions|Gaps|Insertions|Inserts|+ strand %|A minor %|C minor %|G minor %|T minor %|N minor %|Deletion minor %|Insertion minor %
test2|5|A|4|A||1.0|1.87433193885|4|0|0|0|0|0|0|0||0.75|0.25|0.0|0.0|0.0|0.0|0.0|0.0
test2|6|T|4|T/C||1.0|2.94335833285|0|2|0|2|0|0|0|0||0.75|0.0|0.5|0.0|0.0|0.0|0.0|0.0
test2|7|C|4|C|T|1.0|3.45990045591|0|3|0|1|0|0|0|0||0.75|0.0|0.333333333333|0.0|0.0|0.0|0.0|0.0
test2|8|G|4|G||1.0|1.87433193885|0|0|4|0|0|0|0|0||0.75|0.0|0.0|0.25|0.0|0.0|0.0|0.0
'''.replace('|', '\t')

        self.assertEqual(valid, out.getvalue())

    def testBaseCallRegionStrand(self):
        #  confirm that BED strand *isn't* used
        bam = MockBam(['test2'])
        bam.add_read('foo1', 'atcgatcg', '........', 0, 0, cigar='8M')
        bam.add_read('foo2', 'atcgatcg', 'AAAAAAAA', 0, 4, cigar='8M')
        bam.add_read('foo3', 'accgatcg', '########', 0, 4, cigar='8M', is_reverse=True)
        bam.add_read('foo4', 'atcgactgatcg', '############', 0, 0, cigar='12M')

        bed = BedFile(fileobj=StringIO.StringIO('''\
test2|4|7|foo|1|-
'''.replace('|', '\t')))

        #  remember: bed is 0-based, basecall is 1-based, so 4->7, corresponds to 5->8

        out = StringIO.StringIO('')
        ngsutils.bam.basecall.bam_basecall(bam, os.path.join(os.path.dirname(__file__), 'test.fa'), showstrand=True, regions=bed, out=out)

        valid = '''chrom|pos|ref|count|consensus call|minor call|ave mappings|entropy|A|C|G|T|N|Deletions|Gaps|Insertions|Inserts|+ strand %|A minor %|C minor %|G minor %|T minor %|N minor %|Deletion minor %|Insertion minor %
test2|5|A|4|A||1.0|1.87433193885|4|0|0|0|0|0|0|0||0.75|0.25|0.0|0.0|0.0|0.0|0.0|0.0
test2|6|T|4|T/C||1.0|2.94335833285|0|2|0|2|0|0|0|0||0.75|0.0|0.5|0.0|0.0|0.0|0.0|0.0
test2|7|C|4|C|T|1.0|3.45990045591|0|3|0|1|0|0|0|0||0.75|0.0|0.333333333333|0.0|0.0|0.0|0.0|0.0
test2|8|G|4|G||1.0|1.87433193885|0|0|4|0|0|0|0|0||0.75|0.0|0.0|0.25|0.0|0.0|0.0|0.0
'''.replace('|', '\t')

        self.assertEqual(valid, out.getvalue())

    def testBaseCallVariants(self):
        bam = MockBam(['test2'])
        bam.add_read('foo1', 'atcgaccg', '........', 0, 0, cigar='8M')
        bam.add_read('foo2', 'atcgattg', 'AAAAAAAA', 0, 4, cigar='8M')
        bam.add_read('foo3', 'accgattg', '########', 0, 4, cigar='8M')

        # ref: atcgatcgatcg
        #      atcgaccg
        #          atcgattg
        #          accgattg
        #           *    *

        out = StringIO.StringIO('')
        ngsutils.bam.basecall.bam_basecall(bam, os.path.join(os.path.dirname(__file__), 'test.fa'), variants=True, out=out)

        valid = '''chrom|pos|ref|count|consensus call|minor call|ave mappings|entropy|A|C|G|T|N|Deletions|Gaps|Insertions|Inserts
test2|6|T|3|C|T|1.0|2.9097148334|0|2|0|1|0|0|0|0|
test2|11|C|2|T||1.0|4.24102241591|0|0|0|2|0|0|0|0|
'''.replace('|', '\t')

        self.assertEqual(valid, out.getvalue())

    def testBaseCallHetTest(self):
        bam = MockBam(['test2'])
        bam.add_read('foo1', 'atcgaccg', '........', 0, 0, cigar='8M')
        bam.add_read('foo2', 'atcgattg', 'AAAAAAAA', 0, 4, cigar='8M')
        bam.add_read('foo3', 'accgattg', '########', 0, 4, cigar='8M')

        # ref: atcgatcgatcg
        #      atcgaccg
        #          atcgattg
        #          accgattg
        #           *    *

        out = StringIO.StringIO('')
        ngsutils.bam.basecall.bam_basecall(bam, os.path.join(os.path.dirname(__file__), 'test.fa'), variants=True, hettest=True, out=out)

        valid = '''chrom|pos|ref|count|consensus call|minor call|ave mappings|alt. allele freq|entropy|A|C|G|T|N|Deletions|Gaps|Insertions|Inserts
test2|6|T|3|C|T|1.0|0.333333333333|2.9097148334|0|2|0|1|0|0|0|0|
test2|11|C|2|T||1.0|0.0|4.24102241591|0|0|0|2|0|0|0|0|
'''.replace('|', '\t')

        # print '---'
        # print out.getvalue().replace('\t', '|')
        # print '---'
        # print valid.replace('\t', '|')

        self.assertEqual(valid, out.getvalue())

    def testBaseCallMinCount(self):
        bam = MockBam(['test2'])
        bam.add_read('foo1', 'atcgatcg', '........', 0, 0, cigar='8M')
        bam.add_read('foo2', 'atcgatcg', 'AAAAAAAA', 0, 4, cigar='8M')
        bam.add_read('foo3', 'accgatcg', '########', 0, 4, cigar='8M', tags=[('IH', 2), ('HI', 1)])

        out = StringIO.StringIO('')
        ngsutils.bam.basecall.bam_basecall(bam, os.path.join(os.path.dirname(__file__), 'test.fa'), min_count=2, out=out)

        valid = '''chrom|pos|ref|count|consensus call|minor call|ave mappings|entropy|A|C|G|T|N|Deletions|Gaps|Insertions|Inserts
test2|5|A|3|A||1.33333333333|1.65545182665|3|0|0|0|0|0|0|0|
test2|6|T|3|T|C|1.33333333333|2.9097148334|0|1|0|2|0|0|0|0|
test2|7|C|3|C||1.33333333333|4.85048518563|0|3|0|0|0|0|0|0|
test2|8|G|3|G||1.33333333333|1.65545182665|0|0|3|0|0|0|0|0|
test2|9|A|2|A||1.5|1.36179536943|2|0|0|0|0|0|0|0|
test2|10|T|2|T||1.5|4.24102241591|0|0|0|2|0|0|0|0|
test2|11|C|2|C||1.5|4.24102241591|0|2|0|0|0|0|0|0|
test2|12|G|2|G||1.5|1.36179536943|0|0|2|0|0|0|0|0|
'''.replace('|', '\t')
        self.assertEqual(valid, out.getvalue())

    def testHeterzygosity(self):
        self.assertEqual(0.0, ngsutils.bam.basecall._calculate_heterozygosity(10, 5, 5, 0))
        self.assertEqual(0.5, ngsutils.bam.basecall._calculate_heterozygosity(5, 5, 0, 0))
        self.assertEqual(0.1, ngsutils.bam.basecall._calculate_heterozygosity(9, 1, 0, 0))

if __name__ == '__main__':
    unittest.main()
