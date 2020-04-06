#!/bin/bash

set -o pipefail

file=$(find $MD_INPUTDIR -name "*${MD_SUFFIX}" | sed -n ${PBS_ARRAYID}p) #finds the file
name=$(basename ${file%%$MD_SUFFIX}"") #sample ID

echo "Marking Duplicates for sample $name"

java -jar ${PICARD_JAR} MarkDuplicates \
	I=${file} \
    O=$MD_OUTPUTDIR/${name}DupsMarked.bam \
    M=$MD_OUTPUTDIR/${name}DupsMarked_metrics.txt
