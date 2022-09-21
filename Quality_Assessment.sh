#!/bin/bash

set -o pipefail

module load FastQC/0.11.9-Java-11

module load multiqc/1.11-GCCcore-8.3.0-Python-3.8.2

if [[ -d "$QA_INPUTDIR" ]]; then
	for f in `find $QA_INPUTDIR -name "*${SUFFIX}"`; do
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
