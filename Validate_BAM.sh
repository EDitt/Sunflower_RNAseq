#!/bin/bash

set -o pipefail

for file in `find $VB_INPUTDIR -name "*${VB_SUFFIX}"`; do
	name=$(basename ${file%%$VB_SUFFIX}"") #sample ID
	echo "Validating $name BAM file"
	java -jar ${PICARD_JAR} ValidateSamFile \
	I=${file} \
	MODE=SUMMARY \
	O=$VB_OUTPUTDIR/BAMStats #does this append in a loop?
done
