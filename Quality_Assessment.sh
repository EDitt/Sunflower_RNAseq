#!/bin/bash

cd $QA_INPUTDIR

dir=$(find $(pwd -P) -maxdepth 1 -name "*ds*" | sed -n ${PBS_ARRAYID}p)

if [[ -d "$dir" ]]; then
	for f in $dir/*.fastq.gz; do
		if [[ -f "$f" ]]; then
		fastqc -o $QA_OUTPUTDIR \
		-d $QA_TEMP \
		-t 6 $f
	else
		echo "$f is not a valid file"
	fi
done
else
	echo "$dir" is not a valid directory
fi