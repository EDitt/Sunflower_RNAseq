#!/bin/bash

set -o pipefail

if [ "$ANNOTATION_FORMAT" == "GFF3" ]; then
	echo "Using GFF3 annotation to generate a genome index"
	${STAR_FILE} \
	--runThreadN $NTHREAD \
	--runMode genomeGenerate \
	--genomeDir $GEN_DIR \
	--genomeFastaFiles $GEN_FASTA \
	--sjdbGTFtagExonParentTranscript Parent \
	--sjdbGTFfile $GEN_ANN \
	--sjdbOverhang $SPLICE_JUN
elif [ "$ANNOTATION_FORMAT" == "GTF" ]; then
	echo "Using GTF annotation to generate a genome index"
	${STAR_FILE} \
	--runThreadN $NTHREAD \
	--runMode genomeGenerate \
	--genomeDir $GEN_DIR \
	--genomeFastaFiles $GEN_FASTA \
	--sjdbGTFfile $GEN_ANN \
	--sjdbOverhang $SPLICE_JUN
else
	echo "Please specify whether annotation file is in GTF or GFF3 format"
fi