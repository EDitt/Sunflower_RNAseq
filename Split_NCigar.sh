#!/bin/bash

set -o pipefail

file=$(find $SN_INPUTDIR -name "*${SN_SUFFIX}" | sed -n ${PBS_ARRAYID}p) #finds the file
name=$(basename ${file%%$SN_SUFFIX}"") #sample ID

echo "Splitting N-Cigar Reads for sample $name"

gatk SplitNCigarReads \
    -R ${GEN_FASTA} \
    -I ${file} \
    -O $SN_OUTPUTDIR/${name}_split.bam \
    --tmp-dir ${TEMP_DIR}
    