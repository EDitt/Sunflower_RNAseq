#!/bin/bash

set -o pipefail

ID=$(sed -n ${PBS_ARRAYID}p $ID_NAMES) #finds the IDname

if [[ "MB_FILETYPE" == "transcript" ]]; then
	bams=$(find $MB_INPUTDIR -name "$ID[!0-9]*Aligned.toTranscriptome.out.bam") #bam files with the same ID
	echo "Merging $(ls -1 $bams | wc -l) bam files for sample $ID"
	echo "Bam files being merged are: $(ls -1 $bams)"
	samtools merge -nu $MB_OUTPUTDIR/${ID}.merged.bam $bams
elif [[ "MB_FILETYPE" == "genomic" ]]; then
	bams=$(find $MB_INPUTDIR -name "$ID[!0-9]*Aligned.sortedByCoord.out.bam") #bam files with the same ID
	echo "Merging $(ls -1 $bams | wc -l) bam files for sample $ID and then sorting"
	echo "Bam files being merged are: $(ls -1 $bams)"
	samtools merge -u $MB_OUTPUTDIR/${ID}.merged.bam $bams | samtools sort - -o sample.sorted.bam ###finish
else
	echo "Specify whether the bam files are transcript coordinates or genomic coordinates"
	exit 1
fi
