#!/usr/bin/env python

from distutils.core import setup

setup(name='ngsutils',
      version='0.5.1a',
      description='NGSUtils - Various utilities for working with NGS data',
      author='Marcus Breese',
      author_email='mbreese@stanford.edu',
      url='http://ngsutils.org',
      packages=['ngsutils', 'ngsutils.bam', 'ngsutils.bed', 'ngsutils.fastq', 'ngsutils.gtf', 'ngsutils.support', 'ngsutils.ngs'],
#      scripts=['bin/ngsutils', 'bin/fastqutils', 'bin/bamutils', 'bin/bedutils', 'bin/gtfutils'],
      install_requires = ['cython==0.16', 'pysam>=0.4.1', 'coverage>=3.5.3', 'eta>=0.9', 'swalign>=0.2']
     )
