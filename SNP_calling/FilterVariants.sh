 #!/bin/bash

set -o pipefail

 gatk VariantFiltration \
   -R "${GEN_FASTA}" \
   -V input.vcf.gz \
   -O output.vcf.gz \
   --filter-name "my_filter1" \
   --filter-expression "AB < 0.2" \
   --filter-name "my_filter2" \
   --filter-expression "MQ0 > 50"