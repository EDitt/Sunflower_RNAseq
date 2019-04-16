#!/bin/bash

set -o pipefail

for f1 in /scratch/eld72413/Salty_Nut/BamOut/*2068_S*_L001_Aligned.toTranscriptome.out.bam
do
        f2=${f1%%1_Aligned.toTranscriptome.out.bam}"2_Aligned.toTranscriptome.out.bam"
        f3=${f1%%1_Aligned.toTranscriptome.out.bam}"3_Aligned.toTranscriptome.out.bam"
        f4=${f1%%1_Aligned.toTranscriptome.out.bam}"4_Aligned.toTranscriptome.out.bam"
        g1=${f1%%_S*_L001_Aligned.toTranscriptome.out.bam}"-2_S*_L001_Aligned.toTranscriptome.out.bam"
        g2=${g1%%1_Aligned.toTranscriptome.out.bam}"2_Aligned.toTranscriptome.out.bam"
        g3=${g1%%1_Aligned.toTranscriptome.out.bam}"3_Aligned.toTranscriptome.out.bam"
        g4=${g1%%1_Aligned.toTranscriptome.out.bam}"4_Aligned.toTranscriptome.out.bam"
        name=${f1%%-2068_S*_L001_Aligned.toTranscriptome.out.bam}".bam"

samtools merge -nu $name $f1 $f2 $f3 $f4 $g1 $g2 $g3 $g4
done


f1=$(find $RM_INPUTDIR $(pwd -P) -maxdepth 1 -name "*R1_paired.fq.gz" | sed -n ${PBS_ARRAYID}p)