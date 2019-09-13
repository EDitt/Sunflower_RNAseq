#!/bin/bash

set -o pipefail

#cd $AT_INPUTDIR

#dir=$(find $(pwd -P) -maxdepth 1 -name "*ds*" | sed -n ${PBS_ARRAYID}p)

#if [[ -d "$dir" ]]; then
#	name=$(basename ${dir%%-ds.*}"")
#	echo "Trimming $name reads"
#	for f1 in "$dir"/*${FORWARD_NAMING}; do
#		if [[ -f "$f1" ]]; then
#			if [ "$PE" == "True" ]; then
#				f2=${f1%%$FORWARD_NAMING}"$REVERSE_NAMING"
####new
f1=$(find $AT_INPUTDIR -name "*$FORWARD_NAMING" | sed -n ${PBS_ARRAYID}p)
if [[ -f "$f1" ]]; then
	name=$(basename ${f1%%$FORWARD_NAMING}"")
	echo "Trimming $name reads"
	if [ "$PE" == "True" ]; then
		f2=${f1%%$FORWARD_NAMING}"$REVERSE_NAMING"
		java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE \
		-threads 4 \
		-trimlog $AT_OUTPUTDIR/TrimLog_${name} \
		$f1 $f2 \
		$AT_OUTPUTDIR/${name}_R1_paired.fq.gz $AT_OUTPUTDIR/${name}_R1_unpaired.fq.gz \
		$AT_OUTPUTDIR/${name}_R2_paired.fq.gz $AT_OUTPUTDIR/${name}_R2_unpaired.fq.gz \
		ILLUMINACLIP:$ADAPTERFILE:$SEEDMISMATCH:$PALINDROMECLIP:$SIMPLECLIP:$MINADAPTERLEN:$KEEPREADS \
		LEADING:$LEADCUT TRAILING:$TRAILCUT MINLEN:$MINLENGTH
	else
		java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar SE \
		-threads 4 \
		-trimlog $AT_OUTPUTDIR/TrimLog_${name} \
		$f1 \
		$AT_OUTPUTDIR/${name}_R1.fq.gz \
		ILLUMINACLIP:$ADAPTERFILE:$SEEDMISMATCH:$PALINDROMECLIP:$SIMPLECLIP \
		LEADING:$LEADCUT TRAILING:$TRAILCUT MINLEN:$MINLENGTH
	fi
else
	echo "$f1 is not a valid file"
fi
