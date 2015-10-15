#!/bin/bash
RIBO_FILE='tests/data/rat_rpf.bam'
TRANSCRIPTOME_FASTA='tests/data/rat.rna.fna'
RNA_FILE='tests/data/rat_rna.bam'
TRANSCRIPTS=('NM_012875')
COLOR_SCHEMES=('default' 'colorbrewer' 'rgb' 'greyorfs')

for transcript in "${TRANSCRIPTS[@]}"
do
    for scheme in "${COLOR_SCHEMES[@]}"
    do
        # using transcript name for output directory. this might fail!
        riboplot -b $RIBO_FILE -f $TRANSCRIPTOME_FASTA -t $transcript -c $scheme -o ${transcript}_${scheme} -n $RNA_FILE
    done
done

