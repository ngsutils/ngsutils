#!/usr/bin/env python
"""
Extracts exon/intron/junction/promoter information from a refiso file and
writes each annotation to a separate BED format file.

RefIso is a RefFlat format file with a numeric identifier showing which
transcripts are isoforms as the first column.  This can be produced from
KnownGene or RefFlat files.
"""
import sys,os,gzip
from ngs_utils import dictify,gzip_opener

class RefIso(object): 
    def __init__(self,fname):
        self.genes = {}
        
        with gzip_opener(fname) as f:
            for line in f:
                cols = line.rstrip().split('\t')
                iso_id = cols[0]
                if not iso_id in self.genes:
                    self.genes[iso_id] = {}
                    self.genes[iso_id]['_src_'] = {}
                self.genes[iso_id]['_src_'][cols[2]] = dictify(cols[1:],['name','accn','chrom','strand','#txStart','#txEnd','#cdsStart','#cdsEnd','#exonCount','@#exonStarts','@#exonEnds'])
            
        for iso_id in self.genes:
            chrom = None
            name = None
            min_txStart = None
            max_txEnd = None
            strand = None
            chrom = None
            utrs_5 = []
            utrs_3 = []
            exons = []
            introns = []
            promoters = []
            junctions_5 = []
            junctions_3 = []
            for kg_name in self.genes[iso_id]['_src_']:
                src = self.genes[iso_id]['_src_'][kg_name]
                last_end=None
                
                if not strand:
                    strand = src['strand']
                elif strand != src['strand']:
                    strand = '*'
                
                if strand == '+':
                    promoters.append((src['txStart']-2000,src['txStart']))
                elif strand == '-':
                    promoters.append((src['txEnd'],src['txEnd']+2000))


                if not chrom:
                    chrom = src['chrom']
                else:
                    assert chrom == src['chrom']

                if not name:
                    name = src['name']
                else:
                    assert name == src['name']
                
                
                if not min_txStart or src['txStart'] < min_txStart:
                    min_txStart = src['txStart']
                if not max_txEnd or src['txEnd'] > max_txEnd:
                    max_txEnd = src['txEnd']
                
                
                for start,end in zip(src['exonStarts'],src['exonEnds']):
                    if not (start,end) in exons:
                        exons.append((start,end))

                    if last_end is None:
                        if src['cdsStart']!=src['cdsEnd']:
                            utrs_5.append((src['txStart'],start))
                    else:
                        if not (last_end,start) in introns:
                            introns.append((last_end,start))

                        if not (last_end-1,last_end+1) in junctions_5:
                            junctions_5.append((last_end-1,last_end+1))
                        if not (start-1,start+1) in junctions_3:
                            junctions_3.append((start-1,start+1))
                        

                    last_end = end

                if last_end is not None:
                    if src['cdsStart']!=src['cdsEnd']:
                        utrs_3.append((last_end,src['txEnd']))

            self.genes[iso_id]['chrom'] = chrom
            self.genes[iso_id]['name'] = name
            self.genes[iso_id]['exon'] = exons
            self.genes[iso_id]['intron'] = introns
            self.genes[iso_id]['promoter_2k'] = promoters
            if strand == '+' or strand == '*':
                self.genes[iso_id]['utr_5'] = utrs_5
                self.genes[iso_id]['utr_3'] = utrs_3
                self.genes[iso_id]['junction_d'] = junctions_5
                self.genes[iso_id]['junction_a'] = junctions_3
            else:
                self.genes[iso_id]['utr_5'] = utrs_3
                self.genes[iso_id]['utr_3'] = utrs_5
                self.genes[iso_id]['junction_d'] = junctions_3
                self.genes[iso_id]['junction_a'] = junctions_5
                
            self.genes[iso_id]['strand'] = strand
            self.genes[iso_id]['txStart'] = min_txStart
            self.genes[iso_id]['txEnd'] = max_txEnd


    def write_bed(self,out_templ):
        def _write_bed(fobj,chrom,start,end,name,strand):
            fobj.write('%s\t%s\t%s\t%s\t0\t%s\n' % (chrom,start,end,name,strand))
        def _write_bed_ranges(fobj,chrom,name,strand,ranges):
            for start,end in ranges:
                _write_bed(fobj,chrom,start,end,name,strand)
    
        outputs = ['gene','exon','intron','utr_5','utr_3','promoter_2k','junction_d','junction_a']
        ranges = outputs[1:]
        fobjs = {}
        for out in outputs:
            fobjs[out] = open(out_templ % out,'w')
    
        for iso_id in self.genes:
            strand = self.genes[iso_id]['strand']
        
            _write_bed(fobjs['gene'],self.genes[iso_id]['chrom'],self.genes[iso_id]['txStart'],self.genes[iso_id]['txEnd'],self.genes[iso_id]['name'],strand)
        
            for range_name in ranges:
                _write_bed_ranges(fobjs[range_name],self.genes[iso_id]['chrom'],self.genes[iso_id]['name'],strand,self.genes[iso_id][range_name])

        for out in outputs:
            fobjs[out].close()

if __name__ == '__main__':
    out_templ = 'refiso_%s.bed'
    
    if len(sys.argv) not in [2,3] or not os.path.exists(sys.argv[1]):
        print __doc__
        print 'Usage: %s refiso.txt {out_template}\nConverts refiso format annotations into a series of BED files.' % os.path.basename(sys.argv[0])
        sys.exit(1)

    if len(sys.argv) > 2:
        out_templ = '%s_%%s.bed' % sys.argv[2]
        
    sys.stderr.write('Loading refiso...\n')
    refiso = RefIso(sys.argv[1])

    sys.stderr.write('Writing annotations...\n')
    refiso.write_bed(out_templ)
