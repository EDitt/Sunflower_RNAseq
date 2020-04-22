#!/bin/bash

set -o pipefail

gatk --java-options "-Xmx${memory}" \
        HaplotypeCaller \
        -R "${GEN_FASTA}" \
        -I "${sample}" \
        -O "${HC_OUTPUTDIR}/${sample_name}_${current_intvl_name}_RawGLs.g.vcf" \
        --heterozygosity "${HETEROZYGOSITY}" \
        --native-pair-hmm-threads "${num_threads}" \
        --emit-ref-confidence GVCF \
        ${settings}