#!/usr/bin/env python
'''
Tests for bedutils fromvcf
'''

import unittest
import StringIO

import ngsutils.bed.fromvcf


vcftest = StringIO.StringIO('''\
##HEADER
#CHROM~POS~ID~REF~ALT~QUAL~FILTER~INFO~FORMAT~NA00001~NA00002~NA00003
20~14370~rs6054257~G~A~29~PASS~NS=3;DP=14;AF=0.5;DB;H2~GT:GQ:DP:HQ~0|0:48:1:51,51~1|0:48:8:51,51~1/1:43:5:.,.
20~17330~.~TT~AA~3~q10~NS=3;DP=11;AF=0.017~GT:GQ:DP:HQ~0|0:49:3:58,50~0|1:3:5:65,3~0/0:41:3
'''.replace('~', '\t'))


class FromVCFTest(unittest.TestCase):
    def testBasecall(self):
        valid = '''\
20|14369|14370|G/A|29|+
20|17329|17331|TT/AA|3|+
'''.replace('|', '\t')

        out = StringIO.StringIO('')
        ngsutils.bed.fromvcf._bed_fromvcf(vcftest, out=out)
        self.assertEqual(out.getvalue(), valid)


if __name__ == '__main__':
    unittest.main()
