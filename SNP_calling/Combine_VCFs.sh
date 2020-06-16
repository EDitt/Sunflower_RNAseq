 #!/bin/bash

set -o pipefail

mkdir -p "${CV_OUTPUTDIR}/Intermediates"

if [[ -d "$CV_INPUT" ]]; then
	echo "This is a directory"
	find $CV_INPUT -name "*.vcf" | sort -V > ${CV_OUTPUTDIR}/Intermediates/VCF_parts_list.txt
	VCF_list="${CV_OUTPUTDIR}/Intermediates/VCF_parts_list.txt"
elif [[ -f "$CV_INPUT" ]]; then
	echo "This is a file"
	VCF_list="$CV_INPUT"
else
	echo "Please specify a valid directory or list of input VCF files in the config"
fi

# concatenate VCF files
bcftools concat -f "${VCF_list}" > "${CV_OUTPUTDIR}/Intermediates/${PROJECT}_concat.vcf"

# filter out indels
vcftools --vcf "${CV_OUTPUTDIR}/Intermediates/${PROJECT}_concat.vcf" \
--remove-indels \
--recode \
--recode-INFO-all \
--out "${CV_OUTPUTDIR}/Intermediates/${PROJECT}_no_indels"

### can use gatk to count variants
#gatk CountVariants \
#-V "${CV_OUTPUTDIR}/Intermediates/${PROJECT}_concat.vcf"