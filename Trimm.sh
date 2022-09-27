#!/bin/bash

set -o pipefail

module load Trimmomatic/0.39-Java-1.8.0_144 > /dev/null 2>&1

# arrayID (minimum for slurm is zero)
if [[ "${QUEUE}" == "Slurm" ]]; then
	PBS_ARRAYID=$((SLURM_ARRAY_TASK_ID+1))
	echo "Processing array ${SLURM_ARRAY_TASK_ID} corresponding to line/file # ${PBS_ARRAYID}"
else
	echo "Processing array ${PBS_ARRAYID} through PBS queuing system"
fi

if [[ -d "$AT_INPUT" ]]; then #if input is a DIRECTORY
	f1=$(find $AT_INPUT -name "*$FORWARD_NAMING" | sed -n ${PBS_ARRAYID}p) #find the forward reads
elif [[ -f "$AT_INPUT" ]]; then #if input is a FILE
	f1=$(sed -n ${PBS_ARRAYID}p $AT_INPUT)
else
	echo "Please specify a valid directory or list in the config"
fi
if [[ -f "$f1" ]]; then
	name=$(basename ${f1%%$FORWARD_NAMING}"")
	echo "Trimming $name reads"
	if [ "$PE" == "True" ]; then #if data are paired-end
		f2=${f1%%$FORWARD_NAMING}"$REVERSE_NAMING"
		java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE \
		-threads 4 \
		-trimlog $AT_OUTPUTDIR/TrimLog_${name} \
		$f1 $f2 \
		$AT_OUTPUTDIR/${name}R1_paired.fq.gz $AT_OUTPUTDIR/${name}R1_unpaired.fq.gz \
		$AT_OUTPUTDIR/${name}R2_paired.fq.gz $AT_OUTPUTDIR/${name}R2_unpaired.fq.gz \
		ILLUMINACLIP:$ADAPTERFILE:$SEEDMISMATCH:$PALINDROMECLIP:$SIMPLECLIP:$MINADAPTERLEN:$KEEPREADS \
		LEADING:$LEADCUT TRAILING:$TRAILCUT MINLEN:$MINLENGTH
	else #if data are single-end
		java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar SE \
		-threads 4 \
		-trimlog $AT_OUTPUTDIR/TrimLog_${name} \
		$f1 \
		$AT_OUTPUTDIR/${name}R1.fq.gz \
		ILLUMINACLIP:$ADAPTERFILE:$SEEDMISMATCH:$PALINDROMECLIP:$SIMPLECLIP \
		LEADING:$LEADCUT TRAILING:$TRAILCUT MINLEN:$MINLENGTH
	fi
else
	echo "$f1 is not a valid file"
fi
