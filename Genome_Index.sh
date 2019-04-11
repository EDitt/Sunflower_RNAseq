#!/bin/bash

set -o pipefail

/usr/local/apps/eb/STAR/2.6.1c-foss-2016b/bin/STAR \
--runThreadN $NTHREAD \
--runMode genomeGenerate \
--genomeDir $GEN_DIR \
--genomeFastaFiles $GEN_FASTA \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile $GEN_GFF3 \
--sjdbOverhang $SPLICE_JUN
