ngsutils - Utility programs for analysis of next-gen sequencing data
===
Usage
---
Run the `bamutils`, `bedutils`, `fastqutils`, or `sequtils` scripts. These are the main driver scripts 
that setup the appropriate environment (loading pysam) and run the requested command.

For ease of use, it is recommended that you add `bamutils`, `bedutils`, `fastqutils`, and `sequtils` to your $PATH.

Running **`bamutils`** without any arguments will return a list of command available.  Running **`bamutils command`**
will give you the parameters required for that command (or `bedutils`, etc...)

General usage format:

`bamutils command {options} filename`  

`bedutils command {options} filename`  

`fastqutils command {options} filename`  

`sequtils command {options} filename`  

bamutils
---

Scripts for manipulating and analyzing BAM files

DNA

* basecall      - Base caller
* minorallele   - Calls minor alleles

RNA

* cims          - Finds regions of unusual deletions (CLIP-Seq)
* rpkm          - Calculates RPKM/counts for genes/regions/repeats

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


sequtils
---

Scripts for assembling the gene model used in `bedutils` and `bamutils` scripts. The gene model is similar to the UCSC refFlat or KnownGene
tab-delimited format, except that it adds one column to the beginning indicating isoforms. For some organisms this column can be redundant. 
But for others, it is a required step to ensure annotated isoforms are on the same chromosome and overlap. We are calling this format 
**RefIso**. RefIso files can be compiled from UCSC refFlat or KnownGene files. If needed, these can be automatically downloaded for each 
organism.

Building RefIso

* kg           - Given a number of BED files, calculate the number of samples that overlap regions in a reference BED file
* refflat      - Sorts a BED file (in place)

General

* genesize     - Extract the sizes of genes from the RefIso model (genomic and transcript lengths)
* junctions    - Create a library of potential splice-junctions based upon the RefIso model


Installing
---

Checkout the code and run `build_deps.sh`. This will download and compile `pysam` (and if needed, `Cython`). Cython will be installed into your
user-packages folder for Python (usually in `$HOME/.local`). Pysam will be installed in the `ngsutils/ext` folder and loaded as needed by the main
scripts.

If you need read-only access use:
`git clone git://github.iu.edu/mbreese/ngsutils.git`

Requires

* Python 2.6+ (including development packages)
* curl

Recommended

* samtools
* tabix

---

Marcus Breese <mbreese@iupui.edu>  
Center for Computational Biology and Bioinformatics  
Indiana University School of Medicine


&copy;2010-2011 Trustees of Indiana University