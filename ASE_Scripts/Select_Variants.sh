#!/bin/bash

set -o pipefail

gatk \
SelectVariants \
-V $INPUT_VCF \
-O ${OUTPUT_VCF_DIR}/${OUTPUT_VCF_NAME}.vcf \
-sn $PARENT_A \
-sn $PARENT_B \
--exclude-non-variants true

#the -env flag means it won't include non-variant sites
#the -ef flag (now "--exclude-non-variants") means it won't include filtered sites. If this is enabled, sites that have been marked as filtered (i.e. have anything other than `.` or `PASS` in the FILTER field) will be excluded from the output