#!/bin/bash

set -o pipefail

# Get input file (from directory or list)
if [[ -d "$RM_INPUT" ]]; then #if input is a DIRECTORY
	f1=$(find $RM_INPUT $(pwd -P) -maxdepth 1 -name "*$FORWARD" | sed -n ${PBS_ARRAYID}p)
elif [[ -f "$RM_INPUT" ]]; then #if input is a FILE
	f1=$(sed -n ${PBS_ARRAYID}p $RM_INPUT)
else
	echo "Please specify a valid directory or list of files in the config"
fi

name=$(basename ${f1%%$FORWARD}"") #Name of forward sample

### Get reverse read if paired-end, make sure files are valid
if [[ "$PE" == "True" ]]; then
	f2=${f1%%$FORWARD}"$REVERSE"
	if [[ -f $f1 && -f $f2 ]]; then
		echo "Mapping PE reads for sample $name"
	else
		echo "$f1 and $f2 are not both valid files"
		exit 1
	fi
else
	f2=""
	if [[ -f $f1 ]]; then
		echo "Mapping SE reads for sample $name"
	else
		echo "$f1 is not a valid file"
		exit 1
	fi
fi

### Obtain Sample Name + Lane # to add read group information to mapped files
LANE_NUM=$(grep -o "L00[1-4]" <<< $name | cut -c 4) #Obtain Lane #
SAMPLE_NAME=${name%%[!0-9]*}
ID="${SAMPLE_NAME}:${FLOWCELL_NAME}.${LANE_NUM}"
echo "File name indicates the sample name is ${SAMPLE_NAME} and the lane number is ${LANE_NUM}"
echo "The read group ID field will be ${ID}"

###Genomic Coordinate Output
if [[ "$GENOMIC_COORDINATE_BAMSORTED" == "yes" ]]; then #If genomic alignments should be output as sorted BAM files
	FORMAT="BAM SortedByCoordinate"
	echo "Output Genomic Alignments will be sorted BAM files"
else
	FORMAT="SAM"
	echo "Output Genomic Alignments will be unsorted SAM files"
fi

###STAR mapping
if [[ "$RM_PASS" == "first" ]]; then ###first pass mode
	echo "In first pass Mode"
	/usr/local/apps/eb/STAR/2.6.1c-foss-2016b/bin/STAR \
	--runThreadN $RM_NTHREAD \
	--genomeDir $GEN_DIR \
	--readFilesIn $f1 $f2 \
	--readFilesCommand gunzip -c \
	--outFileNamePrefix $CJ_OUTPUTDIR/"$name" \
	--outFilterMismatchNmax $MAX_MIS \
	--outFilterMultimapNmax $MAX_N \
	--outFilterScoreMinOverLread $MINSCORE_READL \
	--outFilterMatchNminOverLread $MINMATCH_READL \
	--outReadsUnmapped $UNMAP_F \
	--outSAMtype SAM \
	--quantMode - \
	--outSAMattrRGline ID:${ID} LB:${SAMPLE_NAME} PL:${PLATFORM} SM:${SAMPLE_NAME} PU:${ID} \
	--outFilterType BySJout \
	--outSJfilterReads Unique ## could change later to be a filtering step
elif [[ "$RM_PASS" == "second" ]]; then ###second pass mode
	echo "In second pass mode using $NUM_JUNCTIONS junction files"
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
	--outSAMattrRGline ID:${ID} LB:${SAMPLE_NAME} PL:${PLATFORM} SM:${SAMPLE_NAME} PU:${ID} \
	--outFilterType BySJout
else
	echo "Please specify in the config file whether this is first or second pass mode"
	exit 1
fi
