#!/bin/bash
# A sample riboplot run with a sorted Ribo-Seq BAM file, FASTA file of 
# the transcriptome and a transcript.

echo "Please run this script from the main package directory else this would fail"
riboplot -b tests/data/5hRPFsorted.bam -f tests/data/zebrafish.fna \
    -t 'gi|148357119|ref|NM_001098396.1|' -n tests/data/5hmRNAsorted.bam

