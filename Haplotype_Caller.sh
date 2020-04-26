#!/bin/bash

set -o pipefail

# Get input bam file(s)
BAMs=$(sed -n ${PBS_ARRAYID}p $HC_INPUT)
name=$(basename ${BAMs%%.*}"")

# Is input a bam file or list of bams?
line_ext=$(basename ${BAMs##*.})

if [[ "$line_ext" == "bam" ]]; then #if this is a list of .bam files
	INPUT="$BAMs"
        echo "Running Haplotype Caller on file $INPUT"
else #if this is a list of lists
	declare -a BAM_LIST
	while read line; do
		BAM_LIST=("${BAM_LIST[@]}" "I=${line}")
	done < $BAMs
        echo "Merging ${#BAM_LIST[@]} bam files for sample $name"
        java -jar ${PICARD_JAR} MergeSamFiles \
        ${BAM_LIST[*]} \
        O=${HC_OUTPUTDIR}/${name}_merged.bam \
        TMP_DIR=${TEMP_DIR}
        INPUT="${HC_OUTPUTDIR}/${name}_merged.bam"
        echo "Running Haplotype caller on merged file $INPUT"
fi

#get number of threads
threads=$(echo "${HC_QSUB}" | grep -oE 'ppn=[[:alnum:]]+' | cut -d '=' -f 2) #code from sequence handling
#get memory
memory="$(echo "${HC_QSUB}" | grep -oE 'mem=[[:digit:]]+' | cut -f 2 -d '=')G" #code from sequence handling

gatk --java-options "-Xmx${memory}" \
        HaplotypeCaller \
        -R "${GEN_FASTA}" \
        -I "${INPUT}" \
        -O "${HC_OUTPUTDIR}/${name}_Raw.g.vcf" \
        -L "${HC_INTERVALS}" \
        --heterozygosity "${HETEROZYGOSITY}" \
        --native-pair-hmm-threads "${threads}" \
        --emit-ref-confidence GVCF \
        --interval-padding 150 \
        --tmp-dir "${TEMP_DIR}"
