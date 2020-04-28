 #!/bin/bash

set -o pipefail

#Interval region
if [ "$GD_SCAFFOLDS" != "false" ] && [ $PBS_ARRAYID -eq $Maxarray ]; then #if scaffolds have been added and this is the last job array
	echo "Processing Scaffolds"
	INTERVAL="${GD_SCAFFOLDS}"
	name="ScaffoldSeq"
	NumScaffolds=$(< $INTERVAL wc -l)
	if [[ "$NumScaffolds" -gt 100 ]]; then
		MergeRule="true"
		echo "More than 100 Scaffolds, using merge-input-intervals flag"
	else
		MergeRule="false"
	fi
else
	INTERVAL=$(grep -v "^@" $IntList | awk '{print $1":"$2"-"$3}' | sed -n ${PBS_ARRAYID}p)
	echo "Processing chromosomal sequence $INTERVAL"
	name="$(echo "${INTERVAL}" | tr -s ':' '_')"
	MergeRule="false"
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
    --tmp-dir "${TEMP_DIR}" \
    --interval-padding 150 \
    --merge-input-intervals "${MergeRule}"
