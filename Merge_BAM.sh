#!/bin/bash

set -o pipefail

ID=$(sed -n ${PBS_ARRAYID}p $ID_NAMES) #finds the IDname

bams=$(find $MB_INPUTDIR -name "$ID[!0-9]*${MB_SUFFIX}") #files with the same ID

echo "Concatenating $(ls -1 $bams | wc -l) bam files for sample $ID"
echo "Bam files being Concatenated are: $(ls -1 $bams)"

#samtools merge -ru $MB_OUTPUTDIR/${ID}.merged.bam $bams
samtools cat -o $MB_OUTPUTDIR/${ID}.merged.bam $bams
