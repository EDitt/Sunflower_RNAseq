#!/bin/bash

set -o pipefail

ID=$(sed -n ${PBS_ARRAYID}p $ID_NAMES) #finds the sample ID

declare -a NEW_BAM_LIST
for file in `find $MB_INPUTDIR -name "$ID[!0-9]*${MB_SUFFIX}"`; do
	NEW_BAM_LIST=("${NEW_BAM_LIST[@]}" "I=${file}")
done

echo "Merging ${#NEW_BAM_LIST[@]} bam files for Sample #${ID}"

java -jar ${PICARD_JAR} MergeSamFiles \
${NEW_BAM_LIST[*]} \
O=${MB_OUTPUTDIR}/${ID}_merged.bam