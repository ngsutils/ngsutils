#!/usr/bin/env python
'''
Tests for bedutils reduce
'''

import unittest
import StringIO

from ngsutils.fastq import FASTQ
import ngsutils.fastq.convertqual


class EncodeTest(unittest.TestCase):
    def testFQRead(self):
        fq = StringIO.StringIO('''\
@foo
ACGTacgtACGT
+
CDEFGHIJKLMN
''')
        out = StringIO.StringIO('')
        ngsutils.fastq.convertqual._fastq_convertqual(FASTQ(fileobj=fq), out=out)

        out.seek(0)
        fqout = FASTQ(fileobj=out)
        read = fqout.fetch().next()
        self.assertEqual(read.name, 'foo')
        self.assertEqual(read.seq, 'ACGTacgtACGT')
        self.assertEqual(read.qual, "$%&'()*+,-./")


if __name__ == '__main__':
    unittest.main()
