#!/bin/bash

set -o pipefail

ID=$(sed -n ${PBS_ARRAYID}p $ID_NAMES) #finds the IDname
bams=$(find $MB_INPUTDIR -name "$ID[!0-9]*") #bam files with the same ID
echo "Merging $(ls -1 $bams | wc -l) bam files for sample $ID"
echo "Bam files being merged are: $(ls -1 $bams)"
samtools merge -nu $MB_OUTPUTDIR/${ID}.merged.bam $bams
#java -jar picard.jar MergeSamFiles \
#I=$list \
#O=$MB_OUTPUTDIR/${ID}.merged.bam
#samtools merge -nu \
#$MB_OUTPUTDIR/${ID}.merged.bam \
#-b $list \