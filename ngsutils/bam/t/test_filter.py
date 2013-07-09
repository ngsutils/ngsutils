#!/usr/bin/env python
'''
Tests for bamutils filter
'''

import os
import unittest

import ngsutils.bam
import ngsutils.bam.filter
from ngsutils.bam.t import MockBam, MockRead


class FilterTest(unittest.TestCase):
    def testUnique(self):
        ''' Unique sequences '''
        read1 = MockRead('foo1', tid=0, pos=1, seq='AAAAAAAAAAT')
        read2 = MockRead('foo2', tid=0, pos=1, seq='AAAAAAAAAAC')
        read3 = MockRead('foo3', tid=0, pos=1, seq='AAAAAAAAAAC')
        read4 = MockRead('foo4', tid=0, pos=1, seq='AAAAAACCCCC')
        read5 = MockRead('foo5', tid=0, pos=10, seq='CAAAAAAAAAA', is_reverse=True)
        read6 = MockRead('foo6', tid=0, pos=10, seq='TAAAAAAAAAA', is_reverse=True)
        read7 = MockRead('foo7', tid=1, pos=10, seq='TAAAAAAAAAA', is_reverse=True)

        uniq = ngsutils.bam.filter.Unique()
        self.assertTrue(uniq.filter(None, read1))
        self.assertTrue(uniq.filter(None, read2))
        self.assertFalse(uniq.filter(None, read2))
        uniq.close()

        uniq2 = ngsutils.bam.filter.Unique(10)
        self.assertTrue(uniq2.filter(None, read1))
        self.assertFalse(uniq2.filter(None, read2))
        self.assertFalse(uniq2.filter(None, read3))
        self.assertTrue(uniq2.filter(None, read4))
        self.assertTrue(uniq2.filter(None, read5))
        self.assertFalse(uniq2.filter(None, read6))
        self.assertTrue(uniq2.filter(None, read7))
        uniq2.close()

    def testUniqueStart(self):
        'Uniq starts'

        read1 = MockRead('foo1', tid=0, pos=1, seq='AAAAAAAAAAT')
        read2 = MockRead('foo2', tid=0, pos=1, seq='AAAAAAAAAAC')
        read3 = MockRead('foo3', tid=0, pos=2, seq='AAAAAAAAAC')

        uniqpos = ngsutils.bam.filter.UniqueStart()

        self.assertTrue(uniqpos.filter(None, read1))
        self.assertFalse(uniqpos.filter(None, read2))
        self.assertTrue(uniqpos.filter(None, read3))

    def testUniqueStartRev(self):
        'Uniq starts Rev'

        read1 = MockRead('foo1', tid=0, pos=1, seq='TAAAAAAAAAA', is_reverse=True, aend=11)
        read2 = MockRead('foo2', tid=0, pos=1, seq='CAAAAAAAAAA', is_reverse=True, aend=11)
        read3 = MockRead('foo3', tid=0, pos=1, seq='AAAAAAAAAC', is_reverse=True, aend=10)
        read4 = MockRead('foo4', tid=0, pos=4, seq='AAAAAAAAAC', is_reverse=True, aend=14)
        read5 = MockRead('foo5', tid=0, pos=4, seq='AAAAAAAAAC', is_reverse=True, aend=15)
        read6 = MockRead('foo6', tid=0, pos=150000, seq='CAAAAAAAAAA', is_reverse=True, aend=150011)
        read7 = MockRead('foo7', tid=0, pos=150000, seq='TAAAAAAAAAA', is_reverse=True, aend=150011)

        uniqpos = ngsutils.bam.filter.UniqueStart()

        self.assertTrue(uniqpos.filter(None, read1))
        self.assertFalse(uniqpos.filter(None, read2))
        self.assertTrue(uniqpos.filter(None, read3))
        self.assertTrue(uniqpos.filter(None, read4))
        self.assertTrue(uniqpos.filter(None, read5))
        self.assertTrue(uniqpos.filter(None, read6))
        self.assertFalse(uniqpos.filter(None, read7))

    def testBlacklist(self):
        'Blacklist'
        tmp_fname = os.path.join(os.path.dirname(__file__), 'tmp_list')
        with open(tmp_fname, 'w') as f:
            f.write('foo1\nfoo2\n')

        read1 = MockRead('foo1')
        read2 = MockRead('foo2')
        read3 = MockRead('foo3')

        blacklist = ngsutils.bam.filter.Blacklist(tmp_fname)

        self.assertFalse(blacklist.filter(None, read1))
        self.assertFalse(blacklist.filter(None, read2))
        self.assertTrue(blacklist.filter(None, read3))

        blacklist.close()
        os.unlink(tmp_fname)

    def testWhitelist(self):
        'Whitelist'
        tmp_fname = os.path.join(os.path.dirname(__file__), 'tmp_list')
        with open(tmp_fname, 'w') as f:
            f.write('foo1\nfoo2\n')

        read1 = MockRead('foo1')
        read2 = MockRead('foo2')
        read3 = MockRead('foo3')

        whitelist = ngsutils.bam.filter.Whitelist(tmp_fname)

        self.assertTrue(whitelist.filter(None, read1))
        self.assertTrue(whitelist.filter(None, read2))
        self.assertFalse(whitelist.filter(None, read3))

        whitelist.close()
        os.unlink(tmp_fname)

    def testIncludeExcludeRegion(self):
        'Include/Exclude region'

        bam = MockBam(['chr1', 'chr2'])

        read1 = MockRead('foo1', tid=0, pos=1, aend=51)
        read2 = MockRead('foo2', tid=0, pos=100, aend=150)
        read3 = MockRead('foo3', tid=0, pos=125, aend=175)
        read4 = MockRead('foo4', tid=0, pos=200, aend=250)
        read5 = MockRead('foo5', tid=1, pos=125, aend=175)

        region = 'chr1:100-150'
        excludefilter = ngsutils.bam.filter.ExcludeRegion(region)

        self.assertTrue(excludefilter.filter(bam, read1))
        self.assertFalse(excludefilter.filter(bam, read2))
        self.assertFalse(excludefilter.filter(bam, read3))
        self.assertTrue(excludefilter.filter(bam, read4))
        self.assertTrue(excludefilter.filter(bam, read5))
        excludefilter.close()

        includefilter = ngsutils.bam.filter.IncludeRegion(region)

        self.assertFalse(includefilter.filter(bam, read1))
        self.assertTrue(includefilter.filter(bam, read2))
        self.assertTrue(includefilter.filter(bam, read3))
        self.assertFalse(includefilter.filter(bam, read4))
        self.assertFalse(includefilter.filter(bam, read5))

        includefilter.close()

    def testIncludeExcludeBED(self):
        'Include/Exclude BED (BED3, BED6)'

        tmp_fname = os.path.join(os.path.dirname(__file__), 'tmp_list')
        with open(tmp_fname, 'w') as f:
            f.write('chr1\t100\t150\nchr1\t200\t250\n')

        bam = MockBam(['chr1', 'chr2'])

        read1 = MockRead('foo1', tid=0, pos=1, aend=51)
        read2 = MockRead('foo2', tid=0, pos=100, aend=150)
        read3 = MockRead('foo3', tid=0, pos=125, aend=175, is_reverse=True)
        read4 = MockRead('foo4', tid=0, pos=200, aend=250, is_reverse=True)
        read5 = MockRead('foo5', tid=1, pos=125, aend=175)

        include3 = ngsutils.bam.filter.IncludeBED(tmp_fname, 'nostrand')

        self.assertFalse(include3.filter(bam, read1))
        self.assertTrue(include3.filter(bam, read2))
        self.assertTrue(include3.filter(bam, read3))
        self.assertTrue(include3.filter(bam, read4))
        self.assertFalse(include3.filter(bam, read5))

        include3.close()

        exclude3 = ngsutils.bam.filter.ExcludeBED(tmp_fname, 'nostrand')

        self.assertTrue(exclude3.filter(bam, read1))
        self.assertFalse(exclude3.filter(bam, read2))
        self.assertFalse(exclude3.filter(bam, read3))
        self.assertFalse(exclude3.filter(bam, read4))
        self.assertTrue(exclude3.filter(bam, read5))

        exclude3.close()

        tmp_fname = os.path.join(os.path.dirname(__file__), 'tmp_list')
        with open(tmp_fname, 'w') as f:
            f.write('chr1\t100\t150\tfoo\t1\t+\nchr1\t200\t250\tfoo\t1\t-\nchr1\t1000\t1250\tfoo\t1\t+\n')

        include6 = ngsutils.bam.filter.IncludeBED(tmp_fname)

        self.assertFalse(include6.filter(bam, read1))
        self.assertTrue(include6.filter(bam, read2))
        self.assertFalse(include6.filter(bam, read3))
        self.assertTrue(include6.filter(bam, read4))
        self.assertFalse(include6.filter(bam, read5))

        include6.close()

        exclude6 = ngsutils.bam.filter.ExcludeBED(tmp_fname)

        self.assertTrue(exclude6.filter(bam, read1))
        self.assertFalse(exclude6.filter(bam, read2))
        self.assertTrue(exclude6.filter(bam, read3))
        self.assertFalse(exclude6.filter(bam, read4))
        self.assertTrue(exclude6.filter(bam, read5))

        exclude6.close()

        os.unlink(tmp_fname)

    def testMismatchRef(self):
        mismatch = ngsutils.bam.filter.MismatchRef(1, os.path.join(os.path.dirname(__file__), 'test.fa'))

        bam = MockBam(['test1', 'chr1'])
        read0 = MockRead('foo0', 'aaaaaaaaaa', tid=1, pos=0, aend=10, cigar='10M')
        read1 = MockRead('foo1', 'aaaaaaaaaa', tid=0, pos=0, aend=10, cigar='10M')
        read2 = MockRead('foo2', 'aaaaaaaaat', tid=0, pos=0, aend=10, cigar='10M')
        read3 = MockRead('foo3', 'aaaaaaaatt', tid=0, pos=0, aend=10, cigar='10M')
        read4 = MockRead('foo4', 'aaaaaaaaaa')

        self.assertRaises(ValueError, mismatch.filter, bam, read0)  # 'chr1' not in FASTA file
        self.assertTrue(mismatch.filter(bam, read1))
        self.assertTrue(mismatch.filter(bam, read2))  # 1 mismatch
        self.assertFalse(mismatch.filter(bam, read3))  # 2 mismatches
        self.assertFalse(mismatch.filter(bam, read4))  # unmapped
        mismatch.close()

    def testMismatch(self):
        mismatch = ngsutils.bam.filter.Mismatch(1)

        bam = MockBam(['test1', ])
        read1 = MockRead('foo1', 'aaaaaaaaaa', tid=0, pos=0, aend=10, cigar='10M', tags=[('NM', 0)],)
        read2 = MockRead('foo2', 'aaaaaaaaat', tid=0, pos=0, aend=10, cigar='10M', tags=[('NM', 1)])
        read3 = MockRead('foo3', 'aaaaaaaatt', tid=0, pos=0, aend=10, cigar='10M', tags=[('NM', 2)])
        read4 = MockRead('foo4', 'aaaaaaaaaa')
        read5 = MockRead('foo5', 'aaaatttaaa', tid=0, pos=0, aend=7, cigar='4M3I3M', tags=[('NM', 3)])
        read6 = MockRead('foo6', 'aaaatttaat', tid=0, pos=0, aend=7, cigar='4M3I3M', tags=[('NM', 4)])

        self.assertTrue(mismatch.filter(bam, read1))
        self.assertTrue(mismatch.filter(bam, read2))  # 1 mismatch
        self.assertFalse(mismatch.filter(bam, read3))  # 2 mismatches
        self.assertFalse(mismatch.filter(bam, read4))  # unmapped

        self.assertTrue(mismatch.filter(bam, read5))  # 1 mismatch (1 indel)
        self.assertFalse(mismatch.filter(bam, read6))  # 2 mismatches (1 indel + 1 mm)

        mismatch.close()

    def testFlag(self):
        read1 = MockRead('foo1', tid=0, pos=0)
        read2 = MockRead('foo2', is_reverse=True, tid=0, pos=0)   # 0x10
        read3 = MockRead('foo3')  # 0x4
        read4 = MockRead('foo4', flag=0x14)

        flag = ngsutils.bam.filter.MaskFlag(0x4)  # fail unmapped

        self.assertTrue(flag.filter(None, read1))
        self.assertTrue(flag.filter(None, read2))
        self.assertFalse(flag.filter(None, read3))
        self.assertFalse(flag.filter(None, read4))
        flag.close()

    def testMapped(self):
        read1 = MockRead('foo1')
        read2 = MockRead('foo2', tid=0, pos=1)

        mapped = ngsutils.bam.filter.Mapped()

        self.assertFalse(mapped.filter(None, read1))
        self.assertTrue(mapped.filter(None, read2))
        mapped.close()

    def testSecondary(self):
        read1 = MockRead('foo1')
        read2 = MockRead('foo2', is_secondary=True)
        read3 = MockRead('foo3', flag=0x104)

        secondary = ngsutils.bam.filter.SecondaryFlag()

        self.assertTrue(secondary.filter(None, read1))
        self.assertFalse(secondary.filter(None, read2))
        self.assertFalse(secondary.filter(None, read3))
        secondary.close()

    def testQCFail(self):
        read1 = MockRead('foo1')
        read2 = MockRead('foo2', is_qcfail=True)
        read3 = MockRead('foo3', flag=0x304)

        qc = ngsutils.bam.filter.QCFailFlag()

        self.assertTrue(qc.filter(None, read1))
        self.assertFalse(qc.filter(None, read2))
        self.assertFalse(qc.filter(None, read3))
        qc.close()

    def testLength(self):
        read1 = MockRead('foo1', 'A')
        read2 = MockRead('foo2', 'AAAAAAAAAA')
        read3 = MockRead('foo3', 'AAAAAAAAAAAAAAAAAAAA')

        minlen = ngsutils.bam.filter.ReadMinLength(10)
        self.assertFalse(minlen.filter(None, read1))
        self.assertTrue(minlen.filter(None, read2))
        self.assertTrue(minlen.filter(None, read3))
        minlen.close()

        maxlen = ngsutils.bam.filter.ReadMaxLength(10)
        self.assertTrue(maxlen.filter(None, read1))
        self.assertTrue(maxlen.filter(None, read2))
        self.assertFalse(maxlen.filter(None, read3))
        maxlen.close()

    def testTags(self):
        read1 = MockRead('foo1', tags=[('ZZ', 1), ('ZY', 'foo')])
        read2 = MockRead('foo2', tags=[('ZZ', 2), ('ZY', 'bar')])
        read3 = MockRead('foo3', tags=[('ZZ', 3), ('ZY', 'baz')])

        lt = ngsutils.bam.filter.TagLessThan('ZZ', 2)
        self.assertTrue(lt.filter(None, read1))
        self.assertFalse(lt.filter(None, read2))
        self.assertFalse(lt.filter(None, read3))
        lt.close()

        lte = ngsutils.bam.filter.TagLessThanEqual('ZZ', 2)
        self.assertTrue(lte.filter(None, read1))
        self.assertTrue(lte.filter(None, read2))
        self.assertFalse(lte.filter(None, read3))
        lte.close()

        gt = ngsutils.bam.filter.TagGreaterThan('ZZ', 2)
        self.assertFalse(gt.filter(None, read1))
        self.assertFalse(gt.filter(None, read2))
        self.assertTrue(gt.filter(None, read3))
        gt.close()

        gte = ngsutils.bam.filter.TagGreaterThanEqual('ZZ', 2)
        self.assertFalse(gte.filter(None, read1))
        self.assertTrue(gte.filter(None, read2))
        self.assertTrue(gte.filter(None, read3))
        gte.close()

        eq = ngsutils.bam.filter.TagEqual('ZZ', 2)
        self.assertFalse(eq.filter(None, read1))
        self.assertTrue(eq.filter(None, read2))
        self.assertFalse(eq.filter(None, read3))
        eq.close()

        eq2 = ngsutils.bam.filter.TagEqual('ZY', 'foo')
        self.assertTrue(eq2.filter(None, read1))
        self.assertFalse(eq2.filter(None, read2))
        self.assertFalse(eq2.filter(None, read3))
        eq2.close()

    def testMismatchDBSNP(self):
        ''' MISSING TEST / EXPERIMENTAL
        TODO: write test
        '''
        pass

    def testMismatchRefDBSNP(self):
        ''' MISSING TEST / EXPERIMENTAL
        TODO: write test
        '''
        pass


if __name__ == '__main__':
    unittest.main()
