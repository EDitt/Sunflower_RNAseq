#!/bin/bash

set -o pipefail

if [[ ! -d "$TEMP_DIR" ]]; then
	echo "Please specify a temporary directory, exiting..."
	exit 1
fi

file=$(find $PB_INPUTDIR -name "*${PB_SUFFIX}" | sed -n ${PBS_ARRAYID}p) #finds the file
name=$(basename ${file%%$PB_SUFFIX}"") #sample ID

echo "Marking Duplicates for sample $name"

java -jar ${PICARD_JAR} MarkDuplicates \
	I="${file}" \
    O="${PB_OUTPUTDIR}/${name}DupsMarked.bam" \
    M="${PB_OUTPUTDIR}/${name}Duplicate_metrics.txt"

#echo "Splitting N-Cigar Reads for sample $name"

#gatk SplitNCigarReads \
#	-R "${GEN_FASTA}" \
#	-I "${TEMP_DIR}/${name}DupsMarked.bam" \
#	-O "${PB_OUTPUTDIR}/${name}processed.bam" \
#	--tmp-dir "${TEMP_DIR}"

#rm $TEMP_DIR/${name}DupsMarked.bam
