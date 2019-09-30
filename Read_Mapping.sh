#!/bin/bash

set -o pipefail

if [[ -d "$RM_INPUTDIR" ]]; then #if input is a DIRECTORY
	f1=$(find $RM_INPUTDIR $(pwd -P) -maxdepth 1 -name "*R1*.fq.gz" | sed -n ${PBS_ARRAYID}p)
elif [[ -f "$RM_INPUTDIR" ]]; then #if input is a FILE
	f1=$(sed -n ${PBS_ARRAYID}p $RM_INPUTDIR)
else
	echo "Please specify a valid directory or list of files in the config"
fi
if [ "$RM_PASS" == "first" ]; then   ###first pass mode
	if [ "$PE" == "True" ]; then
		name=$(basename ${f1%%_R1_paired.fq.gz}"_")
		f2=${f1%%1_paired.fq.gz}"2_paired.fq.gz"
		if [[ -f $f1 && -f $f2 ]]; then
			echo "Mapping PE reads for sample $name in first pass mode"
			/usr/local/apps/eb/STAR/2.6.1c-foss-2016b/bin/STAR \
			--runThreadN $RM_NTHREAD \
			--genomeDir $GEN_DIR \
			--readFilesIn $f1 $f2 \
			--readFilesCommand gunzip -c \
			--outFileNamePrefix $RM_OUTPUTDIR/"$name" \
			--outFilterMismatchNmax $MAX_MIS \
			--outFilterMultimapNmax $MAX_N \
			--outFilterScoreMinOverLread $MINSCORE_READL \
			--outFilterMatchNminOverLread $MINMATCH_READL \
			--outReadsUnmapped $UNMAP_F \
			--quantMode $QUANT
		else
			echo "$f1 and $f2 are not both valid files"
		fi
	elif [ "$PE" == "False" ]; then
		name=$(basename ${f1%%_R1.fq.gz}"_")
		if [[ -f $f1 ]]; then
			echo "Mapping SE reads for sample $name in first pass mode"
			/usr/local/apps/eb/STAR/2.6.1c-foss-2016b/bin/STAR \
			--runThreadN $RM_NTHREAD \
			--genomeDir $GEN_DIR \
			--readFilesIn $f1 \
			--readFilesCommand gunzip -c \
			--outFileNamePrefix $RM_OUTPUTDIR/"$name" \
			--outFilterMismatchNmax $MAX_MIS \
			--outFilterMultimapNmax $MAX_N \
			--outFilterScoreMinOverLread $MINSCORE_READL \
			--outFilterMatchNminOverLread $MINMATCH_READL \
			--outReadsUnmapped $UNMAP_F \
			--quantMode $QUANT
		else
			echo "$f1 is not a valid file"
		fi
	else
		echo "Please specify in the config file whether date is PE (True/False)"
	fi
elif [ "$RM_PASS" == "second" ]; then   ### second pass mode
	if [ "$PE" == "True" ]; then
		name=$(basename ${f1%%_R1_paired.fq.gz}"_")
		f2=${f1%%1_paired.fq.gz}"2_paired.fq.gz"
		if [[ -f $f1 && -f $f2 ]]; then
			echo "Mapping PE reads for sample $name in second pass mode using ${#junctions[@]} junction files"
			/usr/local/apps/eb/STAR/2.6.1c-foss-2016b/bin/STAR \
			--runThreadN $RM_NTHREAD \
			--genomeDir $GEN_DIR \
			--readFilesIn $f1 $f2 \
			--readFilesCommand gunzip -c \
			--outFileNamePrefix $RM_OUTPUTDIR/"$name" \
			--outFilterMismatchNmax $MAX_MIS \
			--outFilterMultimapNmax $MAX_N \
			--outFilterScoreMinOverLread $MINSCORE_READL \
			--outFilterMatchNminOverLread $MINMATCH_READL \
			--outReadsUnmapped $UNMAP_F \
			--quantMode $QUANT \
			--sjdbFileChrStartEnd $JUNCTIONS
		else
			echo "$f1 and $f2 are not both valid files"
		fi
	elif [ "$PE" == "False" ]; then
		name=$(basename ${f1%%_R1.fq.gz}"_")
		if [[ -f $f1 ]]; then
			echo "Mapping SE reads for sample $name in second pass mode using ${#junctions[@]} junction files"
			/usr/local/apps/eb/STAR/2.6.1c-foss-2016b/bin/STAR \
			--runThreadN $RM_NTHREAD \
			--genomeDir $GEN_DIR \
			--readFilesIn $f1 \
			--readFilesCommand gunzip -c \
			--outFileNamePrefix $RM_OUTPUTDIR/"$name" \
			--outFilterMismatchNmax $MAX_MIS \
			--outFilterMultimapNmax $MAX_N \
			--outFilterScoreMinOverLread $MINSCORE_READL \
			--outFilterMatchNminOverLread $MINMATCH_READL \
			--outReadsUnmapped $UNMAP_F \
			--quantMode $QUANT \
			--sjdbFileChrStartEnd $JUNCTIONS
		else
			echo "$f1 is not a valid file"
		fi
	else
		echo "Please specify in the config file whether date is PE (True/False)"
	fi
fi
