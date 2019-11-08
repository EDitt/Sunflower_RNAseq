#!/bin/bash

set -o pipefail

if [ "$PE" == "True" ]; then
	f1=$(find $TQ_INPUTDIR $(pwd -P) -maxdepth 1 -name "*bam" | sed -n ${PBS_ARRAYID}p)
	name=$(basename ${f1%%.bam}"_")
	echo "Calculating expression for sample $name"
	rsem-calculate-expression \
	-p 8 \
	--bam \
	--no-bam-output \
	--paired-end \
	--strandedness $STRAND \
	--forward-prob 0 \
	$f1 \
	$RSEM_ref/$REF_NAME \
	$TQ_OUTPUTDIR/RSEMOut_"$name"
elif [[ "$PE" == "False" ]]; then
	f1=$(find $TQ_INPUTDIR $(pwd -P) -maxdepth 1 -name "*bam" | sed -n ${PBS_ARRAYID}p)
	name=$(basename ${f1%%.bam}"_")
	echo "Calculating expression for sample $name"
	rsem-calculate-expression \
	-p $TQ_NTHREAD \
	--bam \
	--no-bam-output \
	--strandedness $STRAND \
	--forward-prob 0 \
	--fragment-length-mean \
	--fragment-length-sd \
	$f1 \
	$RSEM_ref/$REF_NAME \
	$TQ_OUTPUTDIR/RSEMOut_"$name"
else
	echo "specify whether data are PE or SE"
fi

