#!/bin/bash

set -o pipefail

ID=$(sed -n ${PBS_ARRAYID}p $ID_NAMES) #finds the IDname
bams=$(find $MB_INPUTDIR -name "$ID[!0-9]*") #bam files with the same ID
list=$(ls -1 $bams)
echo "$list"
samtools merge -nu $MB_OUTPUTDIR/${ID}.merged.bam $list
#java -jar picard.jar MergeSamFiles \
#I=$list \
#O=$MB_OUTPUTDIR/${ID}.merged.bam
#samtools merge -nu \
#$MB_OUTPUTDIR/${ID}.merged.bam \
#-b $list \