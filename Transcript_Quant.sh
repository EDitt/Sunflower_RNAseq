#!/bin/bash

set -o pipefail

if [ "$PE" == "True" ]; then
	f1=$(find $TQ_INPUTDIR $(pwd -P) -maxdepth 1 -name "*bam" | sed -n ${PBS_ARRAYID}p)
	name=$(basename ${f1%%.bam}"_")
	echo "Calculating expression for sample $name"
	rsem-calculate-expression \
	-p $TQ_NTHREAD \
	--alignments \
	--no-bam-output \
	--paired-end \
	$f1 \
	$RSEM_ref/$REF_NAME \
	$TQ_OUTPUTDIR/RSEMOut_"$name"
elif [[ "$PE" == "False" ]]; then
	f1=$(find $TQ_INPUTDIR $(pwd -P) -maxdepth 1 -name "*bam" | sed -n ${PBS_ARRAYID}p)
	name=$(basename ${f1%%.bam}"_")
	echo "Calculating expression for sample $name"
	rsem-calculate-expression \
	-p $TQ_NTHREAD \
	--alignments \
	--no-bam-output \
	--fragment-length-mean $FRAG_MEAN_LEN \
	--fragment-length-sd $FRAG_MEAN_SD \
	$f1 \
	$RSEM_ref/$REF_NAME \
	$TQ_OUTPUTDIR/RSEMOut_"$name"
else
	echo "specify whether data are PE or SE"
fi

