#!/bin/bash

set -o pipefail

if [ "$PE" == "True" ]; then
	f1=$(find $TQ_INPUTDIR $(pwd -P) -maxdepth 1 -name "*bam" | sed -n ${PBS_ARRAYID}p)
	name=$(basename ${f1%%.merged.bam}"_")
	echo "Calculating expression for sample $name"
	rsem-calculate-expression \
	-p 8 \
	--bam \
	--no-bam-output \
	--paired-end \
	--forward-prob 0 \
	$f1 \
	$RSEM_ref/$REF_NAME \
	$TQ_OUTPUTDIR/RSEMOut_"$name"
elif [[ "$PE" == "False" ]]; then
	file=$(find $TQ_INPUTDIR $(pwd -P) -maxdepth 1 -name "*bam" | sed -n ${PBS_ARRAYID}p)	
        name=$(basename ${file%%.merged.bam}"_")
        echo "Calculating expression for sample $name"	
	rsem-calculate-expression \
	-p 8 \
	--bam \
	--no-bam-output \
	--single-end \
	$file \
	$RSEM_ref/$REF_NAME \
	$TQ_OUTPUTDIR/RSEMOut_"$name"

fi

