#!/bin/bash

set -o pipefail

if [[ -d "$QA_INPUTDIR" ]]; then
	for f in `find $QA_INPUTDIR -name "*.fastq.gz"`; do
		if [[ -f "$f" ]]; then
			echo "Checking quality of $f"
			fastqc -o $QA_OUTPUTDIR \
			-d $QA_TEMP \
			-t 6 $f
		else
			echo "$f is not a valid file"
		fi
	done
	multiqc $QA_OUTPUTDIR -o $QA_OUTPUTDIR
else
	echo "$QA_INPUTDIR is not a valid directory"
fi
