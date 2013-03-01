#!/usr/bin/env python
'''
Converts a GFF annotation file to GTF, concentrating on the data elements
that are important for gene annotations (gene, exons, CDS, etc). It will also
propertly calculate the transcripts (mRNA) for each gene.

Genes, pseudogenes, mRNA, ncRNA, rRNA, snoRNA, snRNA, and tRNAs are all processed.

The output file by default exports "exon" and "CDS" features. Gene and transcript
features can be exported with the proper arguments.

The input GFF must contain "ID", "Name", and "Parent" attributes.

If there is DNA sequence appended to the end of the file, it must be indicated
with a ##FASTA tag or the sequence must be in FASTA format.

'''

import sys
import os
import ngsutils.support


def usage(msg=''):
    if msg:
        sys.stdout.write('%s\n\n' % msg)

    sys.stdout.write(__doc__)
    sys.stdout.write('''
Usage: gtfutils fromgff {options} filename.gff

Options:
  -rnas            Export RNA features
  -genes           Export gene features
  -err filename    Write lines that couldn't be processed to this file
''')
    sys.exit(1)

valid_features = ['GENE', 'PSEUDOGENE', 'MRNA', 'NCRNA', 'RRNA', 'SNORNA', 'SNRNA', 'TRNA', 'CDS', 'EXON']
valid_rnas = ['MRNA', 'NCRNA', 'SNORNA', 'SNRNA', 'TRNA', 'RRNA']
valid_genes = ['GENE', 'PSEUDOGENE']


class GFFConverter(object):
    def __init__(self, export_gene=False, export_rna=False, out=sys.stdout):
        self.out = out
        self.export_gene = export_gene
        self.export_rna = export_rna

        self.transcript_genes = {}
        self.genenames = {}

        self.queue = []

    def convert_gff(self, fname, error_out=None):
        self.out.write('# converted from %s using NGSUtils\n' % os.path.basename(fname))

        ref = ''
        pos = 0

        done = False

        def callback():
            return '%s:%s, %s' % (ref, pos, len(self.queue))

        def done_callback():
            return done

        for line in ngsutils.support.gzip_reader(fname, callback=callback, done_callback=done_callback):
            if line[:7] == '##FASTA' or line[0] == '>':
                # we are in sequence, so lets just end here...
                done = True
                break
            if line[:2] == '##':
                self.out.write('#%s' % line)
                continue
            elif line[0] == '#':
                self.out.write(line)
                continue

            cols = line.strip().split('\t')
            ref = cols[0]
            pos = cols[3]
            feature = cols[2].upper()

            if feature not in valid_features:
                continue

            self.queue.append(cols)
            self.process_queue()

        self.process_queue()

        if self.queue and error_out:
            with open(error_out, 'w') as f:
                for cols in self.queue:
                    f.write('%s\n' % ('\t'.join(cols)))

    def process_queue(self):
        '''
        We use a queue to do the processing
        '''
        outqueue = []
        # print 'processing queue ##############################'
        for cols in self.queue:
            # print cols
            feature = cols[2].upper()
            attrs = self.get_attrs(cols[8])
            if feature in valid_genes:
                if not self.process_gene(attrs, cols):
                    # print "##### GENE ERROR, adding back to queue"
                    outqueue.append(cols)
            elif feature in valid_rnas:
                if not self.process_rna(attrs, cols):
                    # print "##### RNA ERROR, adding back to queue"
                    outqueue.append(cols)
            else:
                if not self.process_other(attrs, cols):
                    # print "##### EXON/CDS ERROR, adding back to queue"
                    outqueue.append(cols)

        self.queue = outqueue

    def get_attrs(self, attr_col):
        attrs = {}
        for attr in attr_col.split(';'):
            k, v = attr.split('=')
            attrs[k.upper()] = v
        return attrs

    def process_gene(self, attrs, cols):
        geneid = attrs['ID']
        transcriptid = attrs['ID']
        genename = attrs['NAME']
        self.genenames[geneid] = attrs['NAME']

        if self.export_gene:
            self.out_line(geneid, transcriptid, genename, cols)

        return True

    def process_rna(self, attrs, cols):

        geneid = attrs['PARENT']
        if not geneid in self.genenames:
            return False

        genename = self.genenames[geneid]
        transcriptid = attrs['ID']
        self.transcript_genes[transcriptid] = geneid

        if self.export_rna:
            self.out_line(geneid, transcriptid, genename, cols)

        return True

    def process_other(self, attrs, cols):
        transcriptids = attrs['PARENT'].split(',')

        if transcriptids[0] in self.transcript_genes:
            if not transcriptids[0] in self.transcript_genes:
                return False

            geneid = self.transcript_genes[transcriptids[0]]
        else:
            geneid = transcriptids[0]
            transcriptids = [attrs['ID']]

        if not geneid in self.genenames:
            return False

        genename = self.genenames[geneid]

        if len(transcriptids) > 1:
            error = False
            for t in transcriptids:
                if not t in self.transcript_genes:
                    return False

                if self.transcript_genes[t] != geneid:
                    error = True
            if error:
                sys.stderr.write("WARNING: gene_id doesn't match for transcripts: %s\n" % ','.join(transcriptids))

        for transcriptid in transcriptids:
            self.out_line(geneid, transcriptid, genename, cols)

        return True

    def out_line(self, geneid, transcriptid, genename, cols):
        self.out.write('\t'.join(cols[:8]))
        self.out.write('\tgene_id "%s"; transcript_id "%s"; gene_name "%s";\n' % (geneid, transcriptid, genename))


if __name__ == '__main__':
    fname = None

    export_rna = False
    export_gene = False

    errorout = None
    last = None

    for arg in sys.argv[1:]:
        if last == '-err':
            errorout = arg
            last = None
        elif arg in ['-err']:
            last = arg
        elif arg == '-rnas':
            export_rna = True
        elif arg == '-genes':
            export_gene = True
        elif not fname and (os.path.exists(arg) or arg == '-'):
            fname = arg

    if not fname:
        usage('Missing input file, or file not found!')

    GFFConverter(export_gene, export_rna).convert_gff(fname, errorout)
