#!/usr/bin/env python
'''
Filter reads in a FASTQ file. The filtering criteria can be applied as a
batch, allowing you to use more than one criterion at a time. 

'''
import sys,os

from fastq_utils import read_fastq

def fastq_filter(filter_chain,stats_fname=None):
    for name,seq,qual in filter_chain.filter():
        sys.stdout.write("%s\n%s\n+\n%s\n" % (name,seq,qual))
        
    stats = []
    p = filter_chain
    while p:
        stats.insert(0,(p.__class__.__name__,p.kept,p.altered,p.removed))
        p = p.parent

    sys.stderr.write('Criteria\tKept\tAltered\tRemoved\n')
    for name,kept,altered,removed in stats:
        sys.stderr.write('%s\t%s\t%s\t%s\n' % (name,kept,altered,removed))

    if stats_fname:
        with open(stats_fname,'w') as f:
            f.write('Criteria\tKept\tAltered\tRemoved\n')
            for name,kept,altered,removed in stats:
                f.write('%s\t%s\t%s\t%s\n' % (name,kept,altered,removed))
            


class FASTQReader(object):
    def __init__(self,fname):
        self.parent = None
        self.fname = fname
        
        self.altered = 0
        self.removed = 0
        self.kept = 0
        
    def filter(self):
        for tup in read_fastq(fname):
            self.kept += 1
            yield tup

class TrimFilter(object):
    def __init__(self,parent,trim_seq,mismatch_pct,min_filter_len):
        self.parent = parent
        self.trim_seq = trim_seq
        self.mismatch_pct = mismatch_pct
        self.min_filter_len = min_filter_len

        self.altered = 0
        self.removed = 0
        self.kept = 0

    def filter(self):
        for name,seq,qual in self.parent.filter():
            trimmed = False
            
            for i in xrange(len(seq)-self.min_filter_len):
                matches = 0.0
                total = 0
                for s,a in zip(seq[i:],self.trim_seq):
                    total += 1
                    if s == a:
                        matches += 1
            
                if (matches / total) >= self.mismatch_pct:
                    trimmed = True
                    self.altered += 1
                    if len(seq) == len(qual):
                        yield ('%s #trim' % name,seq[:i],qual[:i])
                    else:
                        yield ('%s #trim' % name,seq[:i],qual[:i-1])
                    break
            
            if not trimmed:
                self.kept += 1
                yield ('%s #trim' % name,seq,qual)


class PairedFilter(object):
    def __init__(self,parent):
        self.parent = parent
        self._last = None
        self.altered = 0
        self.removed = 0
        self.kept = 0
    
    def filter(self):
        for tup in self.parent.filter():
            if not self._last:
                self._last = tup
            elif self._last[0] == tup[0]:
                yield self._last
                yield tup
                self._last = None
                self.kept += 2
            else:
                self.removed += 1
                self._last = tup

        if self._last:
            self.removed += 1
            
class QualFilter(object):
    def __init__(self,parent,min_qual,window_size):
        self.parent = parent
        self.min_qual = min_qual
        self.window_size = window_size
        self.altered = 0
        self.removed = 0
        self.kept = 0

    def filter(self):
        for name,seq,qual in self.parent.filter():
            quals = [ord(q)-33 for q in qual] # convert from phred to int[]
            yielded = False
            for i in xrange(len(qual)):
                acc = 0.0
                for q in quals[i:i+self.window_size]:
                    acc += q
                    
                if (acc / self.window_size) < self.min_qual: # truncate here
                    self.altered += 1
                    yielded = True
                    if len(seq) == len(qual): #basespace or colorspace w/o prefix
                        yield('%s #qual' % name, seq[:i+self.window_size-1], qual[:i+self.window_size-1])
                    else:
                        yield('%s #qual' % name, seq[:i+self.window_size], qual[:i+self.window_size-1])
                    break

            if not yielded:
                self.kept += 1
                yield ('%s #qual' % name,seq,qual)

class SuffixQualFilter(object):
    def __init__(self,parent,val):
        self.parent = parent
        self.value = val
        self.altered = 0
        self.removed = 0
        self.kept = 0

    def filter(self):
        for name,seq,qual in self.parent.filter():
            alt = False
            while qual and qual[-1] == self.value:
                alt = True
                qual=qual[:-1]
                seq=seq[:-1]

            if altered:
                self.altered += 1
            else:
                self.kept += 1
                
            yield ('%s #suff' % name,seq,qual)
        
    
class WildcardFilter(object):
    def __init__(self,parent,max_num):
        self.parent = parent
        self.max_num = max_num
        self.altered = 0
        self.removed = 0
        self.kept = 0

    def filter(self):
        for name,seq,qual in self.parent.filter():
            count = 0
            for base in seq:
                if base in '4.N':
                    count += 1

            if count <= self.max_num:
                self.kept += 1
                yield ('%s #wild' % name,seq,qual)
            else:
                self.removed += 1

class SizeFilter(object):
    def __init__(self,parent,min_size):
        self.parent = parent
        self.min_size = min_size
        self.altered = 0
        self.removed = 0
        self.kept = 0

    def filter(self):
        for name,seq,qual in self.parent.filter():
            if len(qual) >= self.min_size:
                self.kept += 1
                yield ('%s #size' % name,seq,qual)
            else:
                self.removed += 1
    

def usage():
    print __doc__
    print """Usage: fastqutils filter {opts} {filters} file.fastq{.gz}
Options:
  -stats filename             Write filter stats out to a file

Filters:
  -wildcard num               Discard reads with too many wildcards (N or .)

  -size minsize               Discard reads that are too short

  -qual minval window_size    Truncate reads (5'->3') where the quality falls
                              below a threshold (floating average over 
                              window_size)

  -suffixqual minval          Trim away bases from the 3' end with low quality 
                              value should be given as a character (Sanger 
                              scale)(a.k.a Illumina B-trim)

  -trim seq mismatch mintrim  Trim away at least [mintrim] bases that match a 
                              [sequence] (3' adaptor) allowing for a certain 
                              [mismatch] percentage (0.0-1.0)
                              
  -paired                     Only keep reads that are correctly paired
"""
    sys.exit(1)
    

if __name__ == '__main__':
    fname = None
    stats_fname = None
    filters_config = []

    last = None
    args = None
    
    for arg in sys.argv[1:]:
        if last == '-wildcard':
            filters_config.append((WildcardFilter,int(arg)))
            last = None
        elif last == '-size':
            filters_config.append((SizeFilter,int(arg)))
            last = None
        elif last == '-suffixqual':
            filters_config.append((SuffixQualFilter,int(arg)))
            last = None
        elif last == '-qual':
            if not args:
                args = [QualFilter,]
                
            if len(args) == 1:
                args.append(int(arg))
            elif len(args) == 2:
                args.append(int(arg))
                filters_config.append(args)
                last = None
                args = None
                
        elif last == '-trim':
            if not args:
                args = [TrimFilter,]

            if len(args) == 1:
                args.append(arg)
            elif len(args) == 2:
                args.append(float(arg))
            elif len(args) == 3:
                args.append(int(arg))   
                filters_config.append(args)
                last = None
                args = None
        elif last == '-stats':
            stats_fname = arg
            last = None
        elif arg in ['-wildcard','-size','-qual','-suffixqual','-trim','-stats']:
            last = arg
        elif arg == '-h':
            usage()
        elif arg == '-paired':
            filters_config.append((PairedFilter,))
        elif not fname and os.path.exists(arg):
            fname = arg
        
    if not fname or not filters_config:
        usage()
        
    chain = FASTQReader(fname)
    for config in filters_config:
        sys.stderr.write('%s\n'% '\t'.join([str(x) for x in config]))
        
        clazz = config[0]
        opts = config[1:]
        chain = clazz(chain,*opts)
    
    fastq_filter(chain)