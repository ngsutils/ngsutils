#!/usr/bin/env python
"""
RefIso is a format which is similar to RefFlat, except that it adds one
additional column.  The first column is a number indicating what genes should
be considered isoforms of each other.

For well annotated organisms, these is somewhat redundant as the gene-name
should be unique among genes.  However, for less well annotated organisms,
or organisms with a lot of gene duplication, this is a required step.
"""

import sys
import os

from ngs_utils import dictify, gzip_opener
from eta import ETA
import pysam


class RefIso(object):
    def __init__(self, filename):
        self._genes = {}
        self._pos = 0
        with gzip_opener(filename) as f:
            for line in f:
                cols = line.rstrip().split('\t')
                if not cols[0] in self._genes:
                    self._genes[cols[0]] = _RefGene(cols[0], cols[1], cols[3], cols[4])
                self._genes[cols[0]].add_transcript(_RefTranscript.from_cols(cols))

    def fsize(self):
        return len(self._genes)

    def tell(self):
        return self._pos

    @property
    def genes(self):
        self._pos = 0
        for gene in self._genes:
            yield self._genes[gene]
            self._pos += 1


class _RefGene(object):
    def __init__(self, iso_id, name, chrom, strand):
        self.iso_id = iso_id
        self.name = name
        self.chrom = chrom
        self.strand = strand
        self.exons = set()
        self._regions = None
        self._transcripts = {}

        self.tx_start = None
        self.tx_end = None

    def add_transcript(self, transcript):
        self._transcripts[transcript.name] = transcript
        if self.tx_start is None or transcript.tx_start < self.tx_start:
            self.tx_start = transcript.tx_start

        if self.tx_end is None or transcript.tx_end > self.tx_end:
            self.tx_end = transcript.tx_end

    @property
    def transcripts(self):
        for t in self._transcripts:
            yield self._transcripts[t]

    @property
    def regions(self):
        if not self._regions:
            starts = []
            ends = []
            accns = []

            for transcript in self._transcripts:
                accns.append(transcript)
                starts.append(self._transcripts[transcript].exon_starts)
                ends.append(self._transcripts[transcript].exon_ends)

            self._regions = calc_regions(self.tx_start,self.tx_end,accns,starts,ends)

        i=0
        for start,end,const,names in self._regions:
            i+=1
            yield (i,start,end,const,names)

class _RefTranscript(object):
    @staticmethod
    def from_cols(cols):
        return _RefTranscript(cols[2],int(cols[5]),int(cols[6]),int(cols[7]),int(cols[8]),[int(x) for x in cols[10].strip(',').split(',')],[int(x) for x in cols[11].strip(',').split(',')])
        
    def __init__(self,name,tx_start,tx_end,cds_start,cds_end,exon_starts,exon_ends):
        self.name = name
        self.tx_start = tx_start
        self.tx_end = tx_end
        self.cds_start = cds_start
        self.cds_end = cds_end
        
        self.exon_starts = exon_starts
        self.exon_ends = exon_ends
        
   
def calc_regions(txStart,txEnd,kg_names,kg_starts,kg_ends):
    '''
        This takes a list of start/end positions (one set per isoform)

        It splits these into regions by assigning each base in the 
        txStart-txEnd range a number.  The number is a bit-mask representing
        each isoform that includes that base in a coding region.  Each
        isoform has a number 2^N, so the bit-mask is just a number.

        Next, it scans the map for regions with the same bit-mask.
        When a different bitmask is found, the previous region (if not intron)
        is added to the list of regions.
        
        Returns a list of tuples:
        (start,end,is_const,names) where names is a comma-separated string 
                                   of accns that make up the region

    '''

    map = [0,] * (txEnd-txStart+1)
    mask = 1

    mask_start_end = {}
    mask_names={}

    for name,starts,ends in zip(kg_names,kg_starts,kg_ends):
        mask_start = None
        mask_end = None
        mask_names[mask] = name
        for start,end in zip(starts,ends):
            if not mask_start:
                mask_start = int(start)
            mask_end = int(end)

            for i in xrange(int(start)-txStart,int(end)-txStart+1):
                map[i] = map[i] | mask

        mask_start_end[mask] = (mask_start,mask_end)
        mask = mask * 2

    last_val = 0
    regions = []
    region_start = 0

    def _add_region(start,end,value):
        rstart = start+txStart
        rend = end+txStart
        const = True
        names = []

        '''
        This code only calls a region alt, if the transcript actually spans this region.
        
        example - two transcripts:
        
        1)        xxxxxx------xxxx------xxx---xxxx--xxxxxxxxx
        2)           xxx------xxxx--xxx-------xxxx--xxxxxxxxxxxxx
        const/alt cccccc------cccc--aaa-aaa---cccc--ccccccccccccc

        I'm not sure if this is a good idea or not...
        
        '''
        
        for mask in mask_start_end:
            mstart,mend = mask_start_end[mask]
            if rstart >= mstart and rend <= mend:
                if value & mask == 0:
                    const = False
                else:
                    names.append(mask_names[mask])

        regions.append((rstart,rend,const,','.join(names)))


        # Alternative code to call everything at ends...
        
        # const = True
        # names = []
        # for mask in mast_start_end:
        #     if value & mask == 0:
        #         const = False
        #     else:
        #         names.append(mask_names[mask])
        # 
        # regions.append((rstart,rend,const,','.join(names)))



    for i in xrange(0,len(map)):
        if map[i]==last_val:
            continue

        if last_val:
            _add_region(region_start,i-1,last_val)

            region_start = i-1
        else:
            region_start = i

        last_val = map[i]

    if last_val:
        _add_region(region_start,i,last_val) # extend by one...

    return regions 

def kg_to_refiso(kgf,xref,iso,extend_window=0):
    iso_kgs = {}
    with gzip_opener(iso) as f:
        for line in f:
            cols = line.rstrip().split('\t')
            if not cols[0] in iso_kgs:
                iso_kgs[cols[0]] = []
            iso_kgs[cols[0]].append(cols[1])

    xref_names = {}
    with gzip_opener(xref) as f:
        for line in f:
            cols = line.rstrip().split('\t')
            xref_names[cols[0]]=cols[4]

    kg_lines = {}
    with gzip_opener(kgf) as f:
        for line in f:
            line=line.rstrip()
            cols = line.split('\t')
            kg_lines[cols[0]]=line

    isos = [ int (x) for x in iso_kgs.keys() ]
    isos.sort()

    iso_count = 0
    for iso in [ str(x) for x in isos ]:
        if len(iso_kgs[iso]) == 1:
            name = xref_names[iso_kgs[iso][0]]
        else:
            names = {}
            for kg in iso_kgs[iso]:
                if not xref_names[kg] in names:
                    names[xref_names[kg]] = 1
                else:
                    names[xref_names[kg]] += 1

            sorted_names = []
            for name in names:
                sorted_names.append((names[name],name))
            sorted_names.sort()
            name = sorted_names[-1][1]

        names = []
        starts = []
        ends = []
        
        cols = {}
        
        for kg in iso_kgs[iso]:
            cols[kg] = kg_lines[kg].split('\t')[:10]
            names.append(kg)
            starts.append(int(cols[kg][3]))
            ends.append(int(cols[kg][4]))

        clusters = split_isoform_clusters(names,starts,ends,extend_window)
        for names,start,end in clusters:
            iso_count += 1
            for kg in names:
                sys.stdout.write('%s\t%s\t' % (iso_count,name))
                sys.stdout.write('\t'.join(cols[kg]))
                sys.stdout.write('\n')
        

def split_isoform_clusters(names,starts,ends,extend_window = 0,verbose=False):
    '''
    given a set of names, start/end positions (txStart/End), this method
    will determine how many isoform clusters there should be.
    
    (Names can be anything really, it isn't assumed to be a string)
    
    It does this by looking for overlapping transcripts.  This
    method assumes all given transcripts are on the same chrom/strand.
    You can also specify an optional wiggle factor for merging transcripts
    together.
    
    Returns a list of tuples ([names,],start,end) that identify the
    isoform clusters
    '''
    clusters = []

    for name,start,end in zip(names,starts,ends):
        clusters.append(([name,],start,end))
    clusters.sort()

    to_remove = True
    altered=True
    r=0
    while altered:
        altered=False
        to_remove = []
        if verbose:
            print r,clusters
            r+=1
            
        for i in xrange(len(clusters)):
            names,start,end = clusters[i]
            if verbose:
                print '->',names
            for j in xrange(len(clusters)):
                if i == j:
                    continue
                jnames,jstart,jend = clusters[j]
                if verbose:
                    print '    ->',jnames,
                
                if (start-extend_window) <= jstart <= (end+extend_window) or (start-extend_window) <= jend <= (end+extend_window):
                    if verbose:
                        print 'merging'
                    # if this start/end is w/in the base start/end, merge the two.
                    names.extend(jnames)
                    new_start = min(start,jstart)
                    new_end = max(end,jend)
                    clusters.append((names,new_start,new_end))
                    to_remove.append(i)
                    to_remove.append(j)
                    break
                elif verbose:
                    print

            if to_remove:
                to_remove.sort()
                for idx in to_remove[::-1]:
                    del clusters[idx]
                to_remove = []
                altered=True
                break

    return clusters


def refflat_to_refiso(fname,extend_window=0):
    genes = {}

    with gzip_opener(fname) as f:
        for line in f:
            cols = line.rstrip().split('\t')
            k = (cols[0],cols[2],cols[3]) # name,chrom,strand
            k2 = (cols[1],cols[2],cols[3],cols[4],cols[5]) # txname,chrom,strand,txstart,txend
            if not k in genes:
                genes[k] = {}
                genes[k]['_src_'] = {}
                genes[k]['_line_'] = {}
                genes[k]['name'] = cols[0]
                genes[k]['chrom'] = cols[2]
                genes[k]['strand'] = cols[3]

            genes[k]['_src_'][k2] = dictify(cols[1:],['name','chrom','strand','#txStart','#txEnd','#cdsStart','#cdsEnd','#exonCount','@#exonStarts','@#exonEnds'])
            genes[k]['_line_'][k2] = line

    iso_count = 0
    for k in genes:
        genename,chrom,strand = k
        names=[]
        starts=[]
        ends=[]

        for k2 in genes[k]['_src_']:
            name,chrom,strand,start,end = k2
            names.append(k2)
            starts.append(int(start))
            ends.append(int(end))
        
        clusters = split_isoform_clusters(names,starts,ends)

        for names,start,end in clusters:
            iso_count += 1
            for name in names:
                sys.stdout.write('%s\t' % (iso_count))
                sys.stdout.write(genes[k]['_line_'][name])

def refiso_junctions(fname,refname,fragment_size=46,min_size=50,out=sys.stdout,max_exons=3):
    '''
    Given a refiso file and a reference genome, it will produce a fasta file 
    representing all possible unique splice junctions within an isoform.
    
    fragement_size - the maximum amount from each side of a splice to include
    
    min_size       - the minimum length of a junction
    
    max_exons      - the maximum number of exons to include in a junction (for small IG exons)
    
    '''
    
    refiso = RefIso(fname)
    ref = pysam.Fastafile(refname)
    
    references = []
    with open('%s.fai' % refname) as f:
        for line in f:
            cols = line.split('\t')
            references.append(cols[0])

    def _extend_junction(seq,name,chrom,exons,counter=1):
        if counter >= max_exons:
            return
        start,end = exons[0]
        frag_end = end
        if end-start > fragment_size:
            frag_end = start+fragment_size

        seq5 = ref.fetch(chrom,start,frag_end)
        newname = '%s,%s-%s' % (name,start,frag_end)
        newseq = seq + seq5
        if len(newseq) >= min_size:
            yield newname,newseq
            return
        elif len(exons) > 1 and counter+1 < max_exons:
            for i in xrange(1,len(exons)):
                for nn_name,nn_seq in _extend_junction(newseq,newname,chrom,exons[i:],counter+1):
                    yield nn_name,nn_seq

    
    eta=ETA(refiso.fsize(),fileobj=refiso)
    junctions = set()
    for gene in refiso.genes:
        if not gene.chrom in references:
            continue
        for txpt in gene.transcripts:
            exons = zip(txpt.exon_starts,txpt.exon_ends)
            # print exons
            if len(exons) > 1000 or gene.name == 'abParts':
                # skip IG hyper / Ab regions
                continue
            for i,(start,end) in enumerate(exons):
                eta.print_status(extra='%s:%s %s #%s' % (gene.chrom,gene.tx_start,gene.name,i))
                if i == len(exons)-1:
                    # con't splice the last exon
                    continue
                frag_start = start
                
                if end-start > fragment_size:
                    frag_start = end-fragment_size
                
                # print '[%s] %s:%s-%s' % (i,gene.chrom,frag_start,end)
                seq3 = ref.fetch(gene.chrom,frag_start,end)
                for j in xrange(len(exons)-i-1):
                    # print '   [%s]' % (j+i+1),
                    # print '%s-%s' % exons[j+i+1]
                    for name,seq in _extend_junction(seq3,'%s:%s-%s' % (gene.chrom,frag_start,end),gene.chrom,exons[j+i+1:]):
                        if not name in junctions:
                            junctions.add(name)
                            out.write('>%s\n%s\n' % (name,seq))
                
    eta.done()

def usage():
    print __doc__
    print """\
Converts KnownGene or RefFlat formatted files to the RefIso format or extracts 
potential splice junction sequences.  All possible splice junctions are 
compiled for each isoform.

Usage:%s kg knownGene.txt{.gz} kgXref.txt{.gz} knownIsoforms.txt{.gz} \
            {extend_window (bp)}
      %s refflat refFlat.txt{.gz} {extend_window (bp)}

      %s junctions refiso.txt{.gz} ref.fasta {-frag size} {-min size}

Options
-------
Conversion:
    extend_window   allow isoform clusters to be merged if they are 
                    extend_window bases away from each other [default 0]

Junctions:
    ref.fasta       reference genome in FASTA/RAGZ format 
                    (samtools indexed ref.fasta.fai req'd)
    -frag size      the number of bases on either side of the [default 46]
                    junction to include
    -min size       the minimum size of a junction            [default 50]
    

""" % (os.path.basename(sys.argv[0]),os.path.basename(sys.argv[0]),os.path.basename(sys.argv[0]))
    sys.exit(1)

if __name__ == '__main__':
    if len(sys.argv) < 2 :
        usage()
    else:
        extend = 0
        if sys.argv[1] == 'kg':
            if len(sys.argv) not in [5,6]:
                usage()
            for arg in sys.argv[2:]:
                if not os.path.exists(arg):
                    usage()
                if len(sys.argv) == 6:
                    extend = int(sys.argv[5])
            kg_to_refiso(sys.argv[2],sys.argv[3],sys.argv[4],extend)

        elif sys.argv[1] == 'refflat':
            if len(sys.argv) not in [3,4]:
                usage()
            if not os.path.exists(sys.argv[2]):
                usage()
            if len(sys.argv) == 4:
                extend = int(sys.argv[3])
            refflat_to_refiso(sys.argv[2],extend)
        elif sys.argv[1] == 'junctions':
            if len(sys.argv) < 4:
                usage()
            if not os.path.exists(sys.argv[2]) or not os.path.exists(sys.argv[3]) or not os.path.exists('%s.fai'%sys.argv[3]):
                usage()
            
            frag_size = 46
            min_size = 50
            
            if len(sys.argv) > 4:
                last = None
                for arg in sys.argv[4:]:
                    if last == '-frag':
                        frag_size = int(arg)
                        last = None
                    elif last == '-min':
                        min_size = int(arg)
                        last = None
                    elif arg in ['-frag','-min']:
                        last = arg
                    else:
                        usage()
            refiso_junctions(sys.argv[2],sys.argv[3],frag_size,min_size)
        else:
            usage()
