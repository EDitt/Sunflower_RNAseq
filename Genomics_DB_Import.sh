 #!/bin/bash

set -o pipefail

#Interval region
if [ "$GD_SCAFFOLDS" != "false" ] && [ $PBS_ARRAYID -eq $Maxarray ]; then #if scaffolds have been added and this is the last job array
	echo "Processing Scaffolds"
	INTERVAL="${GD_SCAFFOLDS}"
	name="ScaffoldSeq"
else
	INTERVAL=$(grep -v "^@" $IntList | awk '{print $1":"$2"-"$3}' | sed -n ${PBS_ARRAYID}p)
	echo "Processing chromosomal sequence $INTERVAL"
	name="$(echo "${INTERVAL}" | tr -s ':' '_')"
fi

#get memory
memory="$(echo "${GD_QSUB}" | grep -oE 'mem=[[:digit:]]+' | cut -f 2 -d '=')"
new_mem_num=$[${memory} - 4]
mem=$(printf "${new_mem_num}g")

gatk --java-options "-Xmx${mem} -Xms${mem}" GenomicsDBImport \
    -R "${GEN_FASTA}" \
	${INPUT} \
	-L "${INTERVAL}" \
    --genomicsdb-workspace-path "${GD_OUTPUTDIR}/gendb_wksp_${name}" \
    --tmp-dir "${TEMP_DIR}"
