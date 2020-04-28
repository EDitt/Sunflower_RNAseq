 #!/bin/bash

set -o pipefail

#Find Interval region
if [ "$GD_SCAFFOLDS" != "false" ] && [ $PBS_ARRAYID -eq $Maxarray ]; then #if scaffolds have been added and this is the last job array
	echo "Genotyping Scaffolds"
	INTERVAL="${GD_SCAFFOLDS}"
	name="ScaffoldSeq"
else
	INTERVAL=$(grep -v "^@" $IntList | awk '{print $1":"$2"-"$3}' | sed -n ${PBS_ARRAYID}p)
	echo "Genotyping chromosomal sequence $INTERVAL"
	name="$(echo "${INTERVAL}" | tr -s ':' '_')"
fi

memory="$(echo "${GV_QSUB}" | grep -oE 'mem=[[:digit:]]+' | cut -f 2 -d '=')g"

cd "${GD_OUTPUTDIR}"

 gatk --java-options "-Xmx${memory}" GenotypeGVCFs \
   -R "${GEN_FASTA}" \
   -L "${INTERVAL}" \
   -V "$gendb_wksp_${name}" \
   -O "${GV_OUTPUTDIR}/${name}.vcf" \
   --tmp-dir "${TEMP_DIR}" \
   --heterozygosity "${HETEROZYGOSITY}" \
   --interval-padding 150
