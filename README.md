ngsutils - Tools for next-gen sequencing analysis
===

Update
---
Development on this tool has largely stopped. For a more up-to-date set of tools that are easier to install and maintain, see https://compgen.io/ngsutilsj


Usage
---
Run the `bamutils`, `bedutils`, `fastqutils`, or `gtfutils` scripts. These are the main driver scripts 
that setup the appropriate environment (loading pysam) and run the requested command.

For ease of use, it is recommended that you add `bamutils`, `bedutils`, `fastqutils`, and `gtfutils` to your $PATH.

Running `bamutils` without any arguments will return a list of commands available.  Running `bamutils command`
will give you the parameters required for that command (or `bedutils`, etc...)

General usage format:

`bamutils command {options} filename`  

`bedutils command {options} filename`  

`fastqutils command {options} filename`  

`gtfutils command {options} filename`  

bamutils
---

Scripts for manipulating and analyzing BAM files

DNA

* basecall      - Base caller
* minorallele   - Calls minor alleles

RNA

* cims          - Finds regions of unusual deletions (CLIP-Seq)
* rpkm          - Calculates RPKM/counts for genes/regions/repeats (Note: gene counts require a *RefIso* file - see `sequtils` below)

General

* convertregion - Converts region mapping to genomic mapping
* expressed     - Finds regions expressed in a BAM file
* extract       - Extract reads from specific regions (BED)
* filter        - Removes reads from a BAM file based on criteria
* merge         - Merge multiple BAM files together
* reads         - Extract the names of reads and their positions
* split         - Splits a BAM file into smaller pieces
* stats         - Calculates simple stats for a BAM file

Conversions

* tobed         - Convert reads to BED6
* tobedgraph    - Convert reads to BedGraph
* tofasta       - Convert reads to FASTA
* tofastq       - Convert reads to FASTQ


bedutils
---

Scripts for manipulating and creating BED files

General

* extend       - Extends BED regions (3' end)
* reduce       - Merges overlapping BED regions
* refcount     - Given a number of BED files, calculate the number of samples that overlap regions in a reference BED file
* sort         - Sorts a BED file (in place)
* stats        - Calculates simple stats for a BED file

Conversions

* fromprimers  - Converts a list of PCR primer pairs to BED regions
* tobedgraph   - BED to BedGraph
* tofasta      - Extract BED regions from a reference FASTA file


fastqutils
---

Scripts for manipulating and creating FASTQ files

General

* convertqual  - Convert a FASTQ file's quality values from Illumina to Sanger scale
* merge        - Merge paired FASTQ files into one file
* names        - Extract the read names from a FASTQ file
* split        - Split one FASTQ file into N number of smaller files
* stats        - Calculates simple stats for a FASTQ file
* trim         - Remove 5' and 3' linkers from each sequence (uses SW alignment)
* truncate     - Truncate the sequence/qual for all reads to a specific length

Conversions

* fromfasta    - Converts FASTA files (with .qual) to FASTQ (basespace or colorspace)
* fromqseq     - Converts Illumina qseq (or export/sorted) files to FASTQ format
* tofasta      - Converts FASTQ to FASTA


gtfutils
---

Scripts for assembling the gene model used in `bedutils` and `bamutils` scripts. The gene model is similar to the UCSC refFlat or KnownGene
tab-delimited format, except that it adds one column to the beginning indicating isoforms. For some organisms this column can be redundant. 
But for others, it is a required step to ensure annotated isoforms are on the same chromosome and overlap. We are calling this format 
**RefIso**. RefIso files can be compiled from UCSC refFlat or KnownGene files. If needed, these can be automatically downloaded for each 
organism. It is possible that this format will be deprecated in the future.

GTF extra annotations

* add_isoform  - Appends isoform annotation from UCSC isoforms file"
* add_reflink  - Appends isoform/name annotation from RefSeq/refLink"
* add_xref     - Appends name annotation from UCSC Xref file"

General

* genesize     - Extract the sizes of genes from the GTF model (genomic and transcript lengths)
* junctions    - Create a library of potential splice-junctions based upon the GTF model

Conversions

* tobed        - Converts a GFF/GTF model to BED format"


Installing
---

Checkout the code and run `make`. This will create a virtualenv folder (env) and install the needed libraries. The only libraries that are
mandatory are *pysam* and *cython*. Cython requires that the Python headers be present on the system. For a linux system this can be
achieved by installing 'python-devel' or similar.

If you need read-only access use:
`git clone git://github.com/ngsutils/ngsutils.git`

Requires

* Python 2.6+ (including development packages)
* virtualenv

Will install

* pysam
* Cython

Recommended

* samtools
* tabix

License
---

    NGSUtils - Tools for next-generation sequencing analysis  
    Copyright (c) 2010-2012 The Trustees of Indiana University
    Copyright (c) 2013-2016 The Board of Trustees of Leland Stanford Junior University

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer listed
      in this license in the documentation and/or other materials
      provided with the distribution.

    - Neither the name of the copyright holders nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.

    The copyright holders provide no reassurances that the source code
    provided does not infringe any patent, copyright, or any other
    intellectual property rights of third parties.  The copyright holders
    disclaim any liability to any recipient for claims brought against
    recipient by any third party for infringement of that parties
    intellectual property rights.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
---

Written and maintained by: Marcus Breese <mbreese@stanford..edu>  
Dept. of Pediatrics, Div. of Hematology/Oncology  
Stanford University School of Medicine


&copy;2010-2012 Trustees of Indiana University  
&copy;2013-2016 The Board of Trustees of Leland Stanford Junior University
