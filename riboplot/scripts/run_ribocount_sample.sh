#!/bin/bash
# A sample ribocount run with a Ribo-Seq sorted BAM format alignment
# and a FASTA file of the transcriptome.
echo "Please run this script from the main package directory else this would fail"
ribocount -b tests/data/5hRPFsorted.bam -f tests/data/zebrafish.fna 

