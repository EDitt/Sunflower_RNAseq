#!/bin/bash

set -o pipefail

# Get sample information
f1=$(find $TQ_INPUTDIR $(pwd -P) -maxdepth 1 -name "*bam" | sed -n ${PBS_ARRAYID}p)
name=$(basename ${f1%%.bam}"_")

# Make sure sample is valid for RSEM
validate=$(rsem-sam-validator ${f1})
if [[ "$validate" =~ "is not valid" ]]; then
	echo "File $name is not valid for RSEM:"
	echo "$validate"
	echo "exiting..."
	exit 1
else
	message=${validate#"........................."}
	echo "RSEM says: $message ...proceeding to quantification"
fi

# RSEM quantification
echo "Calculating expression for sample $name"
if [ "$PE" == "True" ]; then
	rsem-calculate-expression \
	-p $TQ_NTHREAD \
	--alignments \
	--no-bam-output \
	--paired-end \
	--strandedness $STRANDEDNESS \
	--estimate-rspd \
	$f1 \
	$RSEM_ref/$REF_NAME \
	$TQ_OUTPUTDIR/RSEMOut_"$name"
elif [[ "$PE" == "False" ]]; then
	rsem-calculate-expression \
	-p $TQ_NTHREAD \
	--alignments \
	--no-bam-output \
	--strandedness $STRANDEDNESS \
	--estimate-rspd \
	--fragment-length-mean $FRAG_MEAN_LEN \
	--fragment-length-sd $FRAG_MEAN_SD \
	$f1 \
	$RSEM_ref/$REF_NAME \
	$TQ_OUTPUTDIR/RSEMOut_"$name"
else
	echo "specify whether data are PE or SE"
	exit 1
fi

