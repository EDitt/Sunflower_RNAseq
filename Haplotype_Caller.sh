#!/bin/bash

set -o pipefail

# Get input bam file(s)
BAMs=$(sed -n ${PBS_ARRAYID}p $HC_INPUT)
name=$(basename ${BAMs%%.*}"")

# Is input a list of bams or list of lists?
line1=$(sed -n 1p $HC_INPUT)
line1_ext=$(basename ${line1##*.})

if [[ "$line1_ext" == "bam" ]]; then #if this is a list of .bam files
	INPUT="-I $BAMs"
else #if this is a list of lists
	declare -a INPUT
	while read line; do
		INPUT=("${INPUT[@]}" "-I ${line}")
	done < $BAMs
fi

#get number of threads
threads=$(echo "${HC_QSUB}" | grep -oE 'ppn=[[:alnum:]]+' | cut -d '=' -f 2) #code from sequence handling
#get memory
memory="$(echo "${HC_QSUB}" | grep -oE 'mem=[[:digit:]]+' | cut -f 2 -d '=')G" #code from sequence handling

gatk --java-options "-Xmx${memory}" \
        HaplotypeCaller \
        -R "${GEN_FASTA}" \
        ${INPUT} \
        -O "${HC_OUTPUTDIR}/${name}_Raw.g.vcf" \
        -L "${HC_INTERVALS}" \
        --heterozygosity "${HETEROZYGOSITY}" \
        --native-pair-hmm-threads "${threads}" \
        --emit-ref-confidence GVCF \
        --interval-padding 150 \
        --tmp-dir "${TEMP_DIR}"
