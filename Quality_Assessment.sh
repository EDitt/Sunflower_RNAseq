#!/bin/bash

set -o pipefail

dir="${files[${PBS_ARRAYID}]}"
	for f in $dir/*.fastq.gz; do
		fastqc -o $QA_OUTPUTDIR \
		-a $ADAPTERFILE \
		-d $QA_TEMP \
		-t 6 $f
	done