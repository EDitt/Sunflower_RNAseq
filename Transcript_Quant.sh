#!/bin/bash

set -o pipefail

f1=$(find $TQ_INPUTDIR $(pwd -P) -maxdepth 1 -name "*bam" | sed -n ${PBS_ARRAYID}p)
name=$(basename ${f1%%.merged.bam}"_")
rsem-calculate-expression \
-p 8 \
--bam \
--no-bam-output \
--paired-end \
--forward-prob 0 \
$f1 \
$RSEM_ref/$REF_NAME \
$TQ_OUTPUTDIR/RSEMOut_"$name"
