#!/bin/bash

set -o pipefail

# Define Exon Parent Gene and Transcript Tags

if [ "$ANNOTATION_FORMAT" == "GFF3" ]; then
	echo "Using GFF3 annotation to generate a genome index"
	TRANSCRIPT_TAG="Parent"
	if [[ ! -z "$GENE_PARENT" ]]; then
		GENE_TAG="${GENE_PARENT}"
	else
		GENE_TAG="gene_id"
	fi
elif [ "$ANNOTATION_FORMAT" == "GTF" ]; then
	echo "Using GTF annotation to generate a genome index"
	TRANSCRIPT_TAG="transcript_id"
	GENE_TAG="gene_id"
else
	echo "Please specify whether annotation file is in GTF or GFF3 format"
fi
	
${STAR_FILE} \
--runThreadN $NTHREAD \
--runMode genomeGenerate \
--genomeDir $GEN_DIR \
--genomeFastaFiles $GEN_FASTA \
--sjdbGTFtagExonParentTranscript $TRANSCRIPT_TAG \
--sjdbGTFtagExonParentGene $GENE_TAG \
--sjdbGTFfile $GEN_ANN \
--sjdbOverhang $SPLICE_JUN
