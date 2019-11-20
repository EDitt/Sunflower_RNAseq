#!/bin/bash

set -o pipefail

if [[ -d "$RM_INPUT" ]]; then #if input is a DIRECTORY
	f1=$(find $RM_INPUT $(pwd -P) -maxdepth 1 -name "*$FORWARD" | sed -n ${PBS_ARRAYID}p)
elif [[ -f "$RM_INPUT" ]]; then #if input is a FILE
	f1=$(sed -n ${PBS_ARRAYID}p $RM_INPUT)
else
	echo "Please specify a valid directory or list of files in the config"
fi

name=$(basename ${f1%%$FORWARD}"")
f2=${f1%%$FORWARD}"$REVERSE"

###Genomic Coordinate Output
if [[ "$GENOMIC_COORDINATE_BAMSORTED" == "yes" ]]; then #If genomic alignments should be output as sorted BAM files
	FORMAT="BAM SortedByCoordinate"
	echo "Output Genomic Alignments will be sorted BAM files"
else
	FORMAT="SAM"
	echo "Output Genomic Alignments will be unsorted SAM files"
fi

###STAR mapping
if [[ "$RM_PASS" == "first" && "$PE" == "True" ]]; then ###first pass mode, paired-end
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
		--outSAMtype $FORMAT \
		--quantMode $QUANT
	else
		echo "$f1 and $f2 are not both valid files"
		exit 1
	fi
elif [[ "$RM_PASS" == "first" && "$PE" == "False" ]]; then ###first pass mode, single-end
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
		--outSAMtype $FORMAT \
		--quantMode $QUANT
	else
		echo "$f1 is not a valid file"
		exit 1
	fi
elif [[ "$RM_PASS" == "second" && "$PE" == "True" ]]; then ###second pass mode, paired-end
	if [[ -f $f1 && -f $f2 ]]; then
		echo "Mapping PE reads for sample $name in second pass mode using $NUM_JUNCTIONS junction files"
		echo "Junctions are as follows: $JUNCTIONS"
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
		--outSAMtype $FORMAT \
		--quantMode $QUANT \
		--sjdbFileChrStartEnd $JUNCTIONS
	else
		echo "$f1 and $f2 are not both valid files"
		exit 1
	fi
elif [[ "$RM_PASS" == "second" && "$PE" == "False" ]]; then ###second pass mode, single-end
	if [[ -f $f1 ]]; then
		echo "Mapping SE reads for sample $name in second pass mode using $NUM_JUNCTIONS junction files"
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
		--outSAMtype $FORMAT \
		--quantMode $QUANT \
		--sjdbFileChrStartEnd $JUNCTIONS
	else
		echo "$f1 is not a valid file"
		exit 1
	fi
else
	echo "Please specify in the config file whether date is PE (True/False) and whether this is first or second pass mode"
	exit 1
fi
