#!/bin/bash

set -o pipefail

if [ "$ANNOTATION_FORMAT" == "GFF3" ]; then
	rsem-prepare-reference \
	--star  \
	--gff3 $GEN_ANN \
	$GEN_FASTA \
	$RSEM_ref/$REF_NAME
elif [ "$ANNOTATION_FORMAT" == "GTF" ]; then
	rsem-prepare-reference \
	--star  \
	--gtf $GEN_ANN \
	$GEN_FASTA \
	$RSEM_ref/$REF_NAME
else
	echo "Please specify whether annotation file is in GTF or GFF3 format"
fi