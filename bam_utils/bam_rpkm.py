#!/usr/bin/env python
'''
Takes an isoform/gene model, an BAM file and calculates how many reads show 
support of each gene / region.
'''

import sys,os,gzip,re,math
import support.ngs_utils
from support.eta import ETA
from support.refiso import RefIso
import pysam

bam_cigar = ['M','I','D','N','S','H','P']

def calc_rpkm(bam_fname,refiso_name,stranded=True,multiple='complete',normalization='genes',whitelist=None,blacklist=None,coverage=False):
    assert multiple in ['complete','partial','ignore']
    assert normalization in ['genes','total','quartile','none']

    bam = pysam.Samfile(bam_fname,'rb')
    refiso = RefIso(refiso_name)

    genes = []

    all_reads = set()
    eta=ETA(refiso.fsize(),fileobj=refiso)
    for gene in refiso.genes:
        if not gene.chrom in bam.references:
            continue
        eta.print_status(extra=gene.name)

        coding_len = 0
        starts = []
        ends = []

        for num,start,end,const,names in gene.regions:
            coding_len += (end-start)
            starts.append(start)
            ends.append(end)

        count,reads = _fetch_reads(bam,gene.chrom,gene.strand if stranded else None,starts,ends,multiple,False,whitelist,blacklist)
        cols = [gene.name,gene.iso_id,gene.chrom,gene.strand,gene.tx_start,gene.tx_end,coding_len,count]
        if coverage:
            cols.extend(calc_coverage(bam,chrom,strand if stranded else None,[start],[end],whitelist,blacklist))
            eta.print_status(extra='%s coverage' % gene.name)
        genes.append(cols)
        all_reads.update(reads)

    mapped_count = 0

    if normalization == 'genes':
        mapped_count = len(all_reads)
    elif normalization == 'total':
        mapped_count = _find_mapped_count(bam,bam_fname,whitelist)
    elif normalization == 'quartile':
        mapped_count = _find_mapped_count_quartile([gene[-1] for gene in genes])
        
    print "#multiple\t%s" % multiple
    if not stranded:
        print "#nostrand"

    print "#normalization\t%s" % normalization
    headers = "gene isoid chrom strand txStart txEnd coding_len count".split()

    if coverage:
        headers.append('Coverage mean')
        headers.append('Coverage stdev')
        headers.append('Coverage median')
    
    if mapped_count:
        print "#mapped_count\t%s" % mapped_count
        headers.append('RPKM')
        mil_mapped = mapped_count / 1000000.0
    else:
        mil_mapped = 0
    
    print ""
    print '\t'.join(headers)
    
    for gene in genes:
        if mil_mapped:
            if coverage:
                count = gene[-4]
                size = gene[-5]
            else:
                count = gene[-1]
                size = gene[-2]
            gene.append(count / (size/1000.0) / mil_mapped)
        sys.stdout.write('%s\n' % '\t'.join([str(x) for x in gene]))

    eta.done()
    bam.close()

def calc_bed(bam_fname,bed_name,stranded=True,multiple='complete',normalization='genes',whitelist=None,blacklist=None,uniq=False,coverage=False):
    assert multiple in ['complete','partial','ignore']
    assert normalization in ['genes','total','quartile','none']
    
    bam = pysam.Samfile(bam_fname,'rb')
    bed = open(bed_name)

    eta = ETA(os.stat(bed_name).st_size,fileobj=bed)
    all_reads = set()
    regions = []
    
    for line in bed:
        cols = line.strip().split('\t')
        chrom = cols[0]
        start = int(cols[1])
        end = int(cols[2])
        name = cols[3]
        strand = cols[5]

        if not chrom in bam.references:
            continue
        
        eta.print_status(extra=name)
        coding_len = end-start

        count,reads = _fetch_reads(bam,chrom,strand if stranded else None,[start],[end],multiple,False,whitelist,blacklist,uniq)
        cols = [chrom,start,end,0,name,strand,coding_len,count]

        if coverage:
            cols.extend(calc_coverage(bam,chrom,strand if stranded else None,[start],[end],whitelist,blacklist))
            eta.print_status(extra='%s coverage' % name)
        
        regions.append(cols)
        all_reads.update(reads)

    mapped_count = 0
    if normalization == 'genes':
        mapped_count = len(all_reads)
    elif normalization == 'total':
        mapped_count = _find_mapped_count(bam,bam_fname,whitelist,blacklist)
    elif normalization == 'quartile':
        mapped_count = _find_mapped_count_quartile([region[-1] for region in regions])
        

    print "#multiple\t%s" % multiple
    if not stranded:
        print "#nostrand"

    print "#normalization\t%s" % normalization
    headers = "chrom start end name score strand size count".split()

    if coverage:
        headers.append('coverage_mean')
        headers.append('coverage_stdev')
        headers.append('coverage_median')

    if mapped_count:
        print "#mapped_count\t%s" % mapped_count
        headers.append('RPKM')
        mil_mapped = mapped_count / 1000000.0
    else:
        mil_mapped = 0
    
    
    print ""
    print '\t'.join(headers)

    for region in regions:
        if mil_mapped:
            if coverage:
                count = region[-4]
                size = region[-5]
            else:
                count = region[-1]
                size = region[-2]
            region.append(count / (size/1000.0) / mil_mapped)
        print '\t'.join([str(x) for x in region])

    eta.done()
    bam.close()
    bed.close()

def calc_repeat(bam_fname,repeat_fname,stranded=True,multiple='complete',normalization='genes',whitelist=None,blacklist=None):
    assert multiple in ['complete','partial','ignore']
    assert normalization in ['genes','total','none']

    bam = pysam.Samfile(bam_fname,'rb')
    
    with support.ngs_utils.gzip_opener(repeat_fname) as repeat_f:
        repeat_f.next()
        repeat_f.next()
        repeat_f.next()

        eta = ETA(os.stat(repeat_fname).st_size,fileobj=repeat_f)
        all_reads = set()
        repeats = {}
        for line in repeat_f:
            cols = line.strip().split()
            chrom = cols[4]
            start = int(cols[5])-1
            end = int(cols[6])
            strand = '+' if cols[8] == '+' else '-'
            family = cols[10]
            member = cols[9]
            name = '%s|%s' % (family,member)
        
            if not (family,member) in repeats:
                repeats[(family,member)] = {'count':0,'size':0 }

            if not (family,'*') in repeats:
                repeats[(family,'*')] = {'count':0,'size':0 }

            if not chrom in bam.references:
                continue

            eta.print_status(extra=name)
            size = end-start
            repeats[(family,'*')]['size'] += size
            repeats[(family,member)]['size'] += size

            count,reads = _fetch_reads(bam,chrom,strand if stranded else None,[start],[end],multiple,False,whitelist,blacklist)
            repeats[(family,'*')]['count'] += count
            repeats[(family,member)]['count'] += count
            all_reads.update(reads)

    mapped_count = 0
    if normalization == 'genes':
        mapped_count = len(all_reads)
    elif normalization == 'total':
        mapped_count = _find_mapped_count(bam,bam_fname,whitelist,blacklist)


    print "#multiple\t%s" % multiple
    if not stranded:
        print "#nostrand"

    print "#normalization\t%s" % normalization
    headers = "family repeat length count".split()

    if mapped_count:
        print "#mapped_count\t%s" % mapped_count
        headers.append('RPKM')
        mil_mapped = mapped_count / 1000000.0
    else:
        mil_mapped = 0


    print ""
    print '\t'.join(headers)

    keys = []
    for k in repeats:
        keys.append(k)
        
    keys.sort()

    for k in keys:
        if k[1] != '*':
            continue
        cols = [k[0],k[1],repeats[k]['size'],repeats[k]['count']]
        if mil_mapped:
            cols.append(cols[-1] / (cols[-2]/1000.0) / mil_mapped)
        print '\t'.join([str(x) for x in cols])
    
    for k in keys:
        if k[1] == '*':
            continue
        cols = [k[0],k[1],repeats[k]['size'],repeats[k]['count']]
        if mil_mapped:
            cols.append(cols[-1] / (cols[-2]/1000.0) / mil_mapped)
        print '\t'.join([str(x) for x in cols])

    eta.done()
    bam.close()

def calc_alt(bam_fname,refiso_name,stranded=True,multiple='complete',normalization='genes',whitelist=None,blacklist=None):
    assert multiple in ['complete','partial','ignore']
    assert normalization in ['genes','total','none']

    bam = pysam.Samfile(bam_fname,'rb')
    refiso = RefIso(refiso_name)

    lines = []
    all_reads = set()
    all_counts = []
    eta=ETA(refiso.fsize(),fileobj=refiso)
    for gene in refiso.genes:
        if not gene.chrom in bam.references:
            continue
        eta.print_status(extra=gene.name)

        #find const regions
        coding_len = 0
        const_regions = []
        starts = []
        ends = []

        last_const = False
        for num,start,end,const,names in gene.regions:
            starts.append(start)
            ends.append(end)
            coding_len += (end-start)
            if const:
                if not last_const:
                    const_regions.append([])
                const_regions[-1].append((start,end))
                last_const = True
            else:
                last_const = False
        
        # find all reads (rpkm)
        total_count,total_reads = _fetch_reads(bam,gene.chrom,gene.strand if stranded else None,starts,ends,multiple,False,whitelist,blacklist)
        all_reads.update(total_reads)
        all_counts.append(total_count)

        # find const-reads
        const_count = 0
        for const_spans in const_regions:
            starts = []
            ends = []
            for start,end in const_spans:
                starts.append(start)
                ends.append(end)
            count,reads = _fetch_reads(bam,gene.chrom,gene.strand if stranded else None,starts,ends,multiple,True,whitelist,blacklist)
            const_count += count
        
        #find counts for each region
        for num,start,end,const,names in gene.regions:
            count,reads = _fetch_reads(bam,gene.chrom,gene.strand if stranded else None,[start],[end],multiple,False,whitelist,blacklist)
            excl_count,excl_reads = _fetch_reads_excluding(bam,gene.chrom,gene.strand if stranded else None,start,end,multiple,whitelist,blacklist)
            
            # remove reads that exclude this region
            for read in excl_reads:
                if read in reads:
                    reads.remove(read)
                    count = count -1
            
            cols = [gene.iso_id,gene.name,total_count,coding_len,const_count,num,'const' if const else 'alt','',gene.chrom,gene.strand,start,end,end-start,count,excl_count]
            
            other_reads = 0
            for read in total_reads:
                if not read in reads and not read in excl_reads:
                    other_reads += 1
                    
            if other_reads > 0:
                cols.append(float(count-excl_count) / other_reads)
            else:
                cols.append('')
                
            lines.append(cols)

    mapped_count = 0
    if normalization == 'genes':
        mapped_count = len(all_reads)
    elif normalization == 'total':
        mapped_count = _find_mapped_count(bam,bam_fname,whitelist,blacklist)
    elif normalization == 'quartile':
        # since each gene will appear multiple times, we have to store the 
        # counts separately from the cols like above
        mapped_count = _find_mapped_count_quartile(all_counts)

    print "# multiple %s" % multiple
    if not stranded:
        print "# nostrand"

    headers = "iso_id gene total_count total_length".split()
    
    if mapped_count:
        print "#mapped_count\t%s" % mapped_count
        headers.append('RPKM')
        mil_mapped = mapped_count / 1000000.0
    else:
        mil_mapped = 0

    headers.extend("const_count region_num const_alt altEvent chrom strand start end length count excl_count alt_index".split())
    
    print ""
    print '\t'.join(headers)
    
    for cols in lines:
        sys.stdout.write('\t'.join([str(x) for x in cols[:4]]))
        if mil_mapped:
            sys.stdout.write('\t')
            sys.stdout.write(str(cols[2] / (cols[3]/1000.0) / mil_mapped))
            sys.stdout.write('\t')
        sys.stdout.write('%s\n' % '\t'.join([str(x) for x in cols[4:]]))
        

    eta.done()
    bam.close()

def calc_coverage(bam,chrom,strand,starts,ends,whitelist,blacklist):
    coverage = []
    
    for start,end in zip(starts,ends):
        for pileup in bam.pileup(chrom,start,end):
            count = 0
            for pileupread in pileup.pileups:
                if blacklist and pileupread.qname in blacklist:
                    continue
                if not whitelist or pileupread.qname in whitelist:
                    if strand:
                        if strand == '+' and pileupread.alignment.is_reverse:
                            continue
                        elif strand == '-' and not pileupread.alignment.is_reverse:
                            continue
                    
                    if not pileupread.is_del:
                        count += 1

            coverage.append(count)
    
    mean,stdev = mean_stdev(coverage)
    coverage.sort()
    if len(coverage) % 2 == 1:
        median = coverage[len(coverage)/2]
    else:
        a = coverage[(len(coverage)/2)-1]
        b = coverage[(len(coverage)/2)]
        median = (a + b) / 2
    return mean,stdev,median
def mean_stdev(l):
    acc = 0
    for el in l:
        acc += el
    
    mean = float(acc) / len(l)
    acc = 0
    for el in l:
        acc += (el-mean)**2
    
    stdev = math.sqrt(float(acc) / (len(l)-1))
    
    return (mean,stdev)


def _find_mapped_count_quartile(counts,min_pct=0.0,max_pct=0.75):
    counts.sort()
    filtered = [ x for x in counts if x > 0 ]
    acc = 0
    min_i = int(min_pct * len(filtered))
    max_i = int(max_pct * len(filtered))
    
    for count in filtered[min_i:max_i]:
        acc += count

    return acc
    

def _find_mapped_count(bam,bam_fname,whitelist=None,blacklist=None):
    sys.stderr.write('Finding number of mapped reads\n')
    bam.seek(0)
    eta = ETA(0,bamfile=bam)
    mapped_count=0
    names=set()
    for read in bam.fetch():
        if blacklist and read.qname in blacklist:
            continue
        if not whitelist or read.qname in whitelist:
            eta.print_status(extra=mapped_count,bam_pos=(read.rname,read.pos))
            if not read.is_unmapped and not read.qname in names:
                mapped_count+=1
                names.add(read.qname)
    eta.done()
    bam.seek(0)
    sys.stderr.write("%s mapped reads\n" % mapped_count)
    return mapped_count
    
def _calc_read_regions(read):
    regions = []
    start = read.pos
    end = read.pos
    for op,length in read.cigar:
        if op == 0:
            end += length
        elif op == 1:
            pass
        elif op == 2:
            pos += length
        elif op == 3:
            regions.append((start,end))
            end += length
            start = end
        
    regions.append((start,end))
    
    return regions
            
    
    
def _fetch_reads_excluding(bam,chrom,strand,start,end,multiple,whitelist=None,blacklist=None):
    '''
    Find reads that exclude this region.  
    
    Example:  This read excludes region 'B'
    
    ------AAAAAAAAA---------BBBBBBBBB---------------------CCCCCCCC---------
              +++++.......................................+++
    
    '''
    
    reads = set()
    count = 0
    
    for read in bam.fetch(chrom,start,end):
        if not strand or (strand=='+' and not read.is_reverse) or (strand=='-' and read.is_reverse):
            excl = True
            for s,e in _calc_read_regions(read):
                if start <= s <=end or start <= e <= end:
                    excl = False
                    break
            if excl:
                reads.add(read.qname)
                count += 1
    return count,reads
            
    
def _fetch_reads(bam,chrom,strand,starts,ends,multiple,exclusive,whitelist=None,blacklist=None,uniq=False):
    '''
    Find reads that match within the given regions...
    
    starts and ends are lists of start and end positions for the various regions
    we are interested in.
    
    Multiple dictates how reads that map to multiple locations are handled.  
    Possible values:
        'complete' - a full 1.0 is added to the count, regardless of how many other places
                     this read maps.
        'partial'  - a fraction of a point is added to the count, based upon how many
                     other places this read maps.  For example, if a read mapped to 3 places
                     on the genome, it would add only 1/3 to the read count.
        'ignore'   - reads that map to multiple places on the genome are ignored.
                     
    Exclusive is a boolean.  If true, then reads must start AND end w/in the given regions in order to
    count.  If false, the read must just be present in the region (start or end).  This is helpful
    in finding reads that map to constant regions as opposed to alternatively expressed regions.
                     
    Whitelist is used to restrict the counting to only reads in the whitelist. 
    (benchmarking only)
    
    uniq will only include reads with uniq starting positions (strand specific)
    
    '''
    if not chrom in bam.references:
        return 0
    
    reads = set()
    start_pos = set()
    count=0
    
    for s,e in zip(starts,ends):
        for read in bam.fetch(chrom,s,e):
            if blacklist and read.qname in blacklist:
                continue
            if not whitelist or read.qname in whitelist:
                if read.is_reverse:
                    k=(read.aend,'-')
                else:
                    k=(read.pos,'+')

                if uniq and k in start_pos:
                        continue
                
                if not strand or (strand=='+' and not read.is_reverse) or (strand=='-' and read.is_reverse):
                    if exclusive:
                        start_ok=False
                        end_ok=False
                    
                        for s1,e1 in zip(starts,ends):
                            if s1 <= read.pos <= e1:
                                start_ok = True
                                break

                        for s1,e1 in zip(starts,ends):
                            if s1 <= read.aend <= e1:
                                end_ok = True
                                break
                        if not start_ok or not end_ok:
                            continue
                        
                    if not read.qname in reads:
                        start_pos.add(k)
                        reads.add(read.qname)
                        if read.tags and 'IH' in read.tags:
                            ih=int(read.opt('IH'))
                        else:
                            ih = 1
                        if ih==1 or multiple=='complete':
                            count += 1
                        elif multiple == 'partial':
                            count += (1.0/ih)
                        else:
                            pass #ignore...
    return count,reads

def usage():
    print __doc__
    print """\
Usage: %s [gene|alt|bed|repeat] annotation_file {opts} bamfile

Common options:
    -nostrand          ignore strand in counting reads
    -norm=<value>      how to normalize RPKM calc
    -coverage          calculate average coverage for all exons
    -multiple=<value>  how to handle reads that map to multiple locations
    -whitelist=<file>  file containing a white-list of read names
                       only these read-names will be used in the calcs
    -blacklist=<file>  file containing a black-list of read names
                       these read-names will not be used in the calcs

[gene]
    Calculate the number of reads that map within the coding regions of each 
    gene. If [matches] is given, an RPKM calculation is also performed, 
    yielding the normalized RPKM value for each gene.

    Annotation: RefIso
    Calculates: # reads, RPKM, coverage


[alt]
    Calculate the number of reads that map to each expressed region
    for all genes. Also, for each gene, the reads mapping to consecutively 
    constant regions are also found.  With these two numbers an alternative
    index is calculated (# reads in a region / # consec. const. reads).
    
    Regions can be exons, or parts of exons, depending on splicing (determined
    by RefIso annotation file).

    Annotation: RefIso
    Calculates: # reads, RPKM, # const region reads, # region reads, alt-index
    
    Alt-index =     (# region reads) - (# region excluding reads)
                 ------------------------------------------------
                            (# NOT region reads)


[bed]
    Calculates the number of reads in a region, where the region is defined 
    in a BED6 formated file: chrom, start (0), end, name, score, strand.

    bed opts:
       -uniq               only count unique starting positions

    Annotation: BED file
    Calculates: # reads, RPKM, coverage


[repeat]
    Calculates the number of reads that map to various repeat regions in the 
    genome.  Repeat regions are defined by repeatmasker.  Output is the number
    of reads that map to each class and type of repeat.

    Annotation: RepeatMasker output file
    Calculates: # reads, RPKM


Possible values for {-norm}:
    genes    - count all reads that map to a gene/region
               (default)
    total    - count all reads that map anywhere on
               the genome
    quartile - count all the reads in the lower 75%%
               of all genes/regions (not available for repeats)
    none     - don't calc RPKM

Possible values for {-multiple}:
    complete - add to the counts of all genes
               (default)
    ignore   - don't add to the count of any genes
    partial  - add a fractional cound to all genes
               (1/number of matches, ex: 1/3)

""" % (os.path.basename(sys.argv[0]),)

if __name__ == '__main__':
    try:
        opts,(cmd,annotation,bam) = support.ngs_utils.parse_args(sys.argv[1:],{'nostrand':False,'coverage':False,'multiple':'complete','norm':'genes','whitelist':None,'blacklist':None,'uniq':False})
    except:
        usage()
        sys.exit(-1)
    
    if not cmd in ['gene','alt','bed','repeat'] or not bam or not annotation or not os.path.exists(bam) or not os.path.exists(annotation):
        usage()
        sys.exit(-1)

    whitelist = None
    if opts['whitelist'] and os.path.exists(opts['whitelist']):
        whitelist = []
        f = open(opts['whitelist'])
        for line in f:
            whitelist.append(line.strip())
        f.close()

    blacklist = None
    if opts['blacklist'] and os.path.exists(opts['blacklist']):
        whitelist = []
        f = open(opts['blacklist'])
        for line in f:
            blacklist.append(line.strip())
        f.close()

    
    if cmd == 'gene':
        calc_rpkm(bam,annotation,not opts['nostrand'],opts['multiple'],opts['norm'],whitelist,blacklist,opts['coverage'])
    elif cmd == 'bed':
        calc_bed(bam,annotation,not opts['nostrand'],opts['multiple'],opts['norm'],whitelist,blacklist,opts['uniq'],opts['coverage'])
    elif cmd == 'alt':
        calc_alt(bam,annotation,not opts['nostrand'],opts['multiple'],opts['norm'],whitelist,blacklist)
    elif cmd == 'repeat':
        calc_repeat(bam,annotation,not opts['nostrand'],opts['multiple'],opts['norm'],whitelist,blacklist)

