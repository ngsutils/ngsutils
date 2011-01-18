ngsutils - Utility programs for analysis of next-gen sequencing data
===

bamutils
---

Commands available with bamutils:

* basecall      - Base caller
* cims_finder   - Finds regions of unusual deletions (CLIP-Seq)
* expressed     - Finds regions expressed in a BAM file
* filter        - Removes reads from a BAM file based on criteria
* region\_to_ref - Converts region mapping to genomic mapping
* rpkm          - Calculates RPKM/counts for genes/regions/repeats
* split         - Splits a BAM file into smaller pieces
* stats         - Calculates simple stats for a BAM file


bedutils
---

Commands available with bedutils:

* extend    - Extends BED regions (3' end)
* reduce    - Merges overlapping BED regions
* ref_count - Given a number of BED files, calculate the number of samples that overlap regions in a reference BED file
* sort      - Sorts a BED file (in place)
