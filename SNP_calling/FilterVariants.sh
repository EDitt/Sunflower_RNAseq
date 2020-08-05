 #!/bin/bash

set -o pipefail

mkdir -p "${CV_OUTPUTDIR}/Intermediates"

gatk IndexFeatureFile \
     -F "${VCF_INPUT}"

gatk VariantFiltration \
	-R "${GEN_FASTA}" \
   	-V "${VCF_INPUT}" \
   	--cluster-window-size 35 \
   	--cluster-size 3 \
   	-filter-name FS --filter-expression "FS > ${FSV}" \
   	-filter-name QD --filter-expression "QD < ${QDV}" \
   	-O "${CV_OUTPUTDIR}/Intermediates/${PROJECT}_filters.vcf"

vcftools --vcf "${CV_OUTPUTDIR}/Intermediates/${PROJECT}_filters.vcf" \
--remove-filtered-all \
--recode \
--recode-INFO-all \
--out "${CV_OUTPUTDIR}/${PROJECT}_filtered"


