#!/bin/bash

set -o pipefail

#   Where is 'sunflower_rnaseq' located?
SUNFLOWER_RNASEQ=$(pwd -P)

ROUTINE="$1" # What routine are we running?
CONFIG="$2" # Where is our config file?

   #If the specified config exists
if [[ -f "${CONFIG}" ]]
then
    source "${CONFIG}" # Source it, providing parameters and software
else # If it doesn't
    echo "Please specify a valid config file." >&2 # Print error message
    exit 1 # Exit with non-zero exit status
fi

#   Where do we output the standard error and standard output files?
ERROR="${SUNFLOWER_RNASEQ}"/ErrorFiles/"${PROJECT}"
mkdir -p "${ERROR}"


#   Run Sunflower_RNAseq
case "${ROUTINE}" in
    1 | Quality_Assessment)
        echo "$(basename $0): Assessing quality..." >&2
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Quality_Assessment.sh" | qsub -l "${QA_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Quality_Assessment
        ;;
    2 | Adapter_Trimming)
        echo "$(basename $0): Trimming Adapters..." >&2
        if [[ -d "$AT_INPUT" ]]; then #if input is a directory
            echo "$AT_INPUT is a directory"
            declare -a files #an array of files
            for f1 in `find $AT_INPUT -name "*$FORWARD_NAMING"`; do
                if [[ -f "$f1" ]]; then
                    files=("${files[@]}" "$f")
                else
                    echo "Please specify a path to valid files in the config file"
                fi
            done
            Maxarray=${#files[@]}
        elif [[ -f "$AT_INPUT" ]]; then #if input is a file
            echo "$AT_INPUT is a file"
            Maxarray=$(< $AT_INPUT wc -l)
        else
            echo "Please specify a valid directory or list in the config"
        fi
        echo "Max array index is ${Maxarray}">&2
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Trimm.sh" | qsub -l "${AT_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Adapter_Trimming -t 1-"${Maxarray}"
        ;;
    3 | Genome_Index)
        if [[ -f "$JUNCTIONS" ]]; then #if regenerating genome index from novel junctions
            echo "$(basename $0): Regenerating genome index with novel junctions discovered in first mapping pass"
        else
            echo "$(basename $0): Generating a genome index from ${ANNOTATION_FORMAT} annotation..." >&2
            echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Genome_Index.sh" | qsub -l "${GI_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Genome_Index
        fi
        ;;
    4 | Collect_Junctions)
        if [[ ! -d "$GEN_DIR" ]]; then #if genome directory is not specified
            echo "Please specify a valid filepath to the genome directory, exiting..."
            exit 1
        fi
        export RM_PASS="first"
        declare -a files #an array of files
        if [[ -d "$RM_INPUT" ]]; then #if input is a directory
            echo "$RM_INPUT is a directory"
            for f in `find $RM_INPUT -name "*$FORWARD"`; do
                if [[ -f "$f" ]]; then
                    files=("${files[@]}" "$f")
                else
                    echo "$f is not a file"
                fi
            done
            Maxarray=${#files[@]}
        elif [[ -f "$RM_INPUT" ]]; then #if input is a file
            echo "$RM_INPUT is a file"
            Maxarray=$(< $RM_INPUT wc -l)
        else
            echo "Please specify a valid directory or list in the config"
            exit 1
        fi
        if [ "$PE" == "True" ]; then
            echo "$(basename $0): Running 1st-pass Mapping to Identify Novel Splice Junctions using PE Reads..." >&2
        elif [ "$PE" == "False" ]; then
            echo "$(basename $0): Running 1st-pass Mapping to Identify Novel Splice Junctions using SE Reads..." >&2
        else
            echo "Please specify in the config file whether data is PE (True/False), exiting..."
            exit 1
        fi
        echo "Max array index is ${Maxarray}">&2
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Read_Mapping.sh" | qsub -l "${RM_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Collect_Junctions  -V -t 1-"${Maxarray}"
        ;;
    5 | Filter_Junctions)
        echo "$(basename $0): Filtering and concatenating junctions for mapping..." >&2
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Filter_Junctions.sh" | qsub -l "${FJ_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Filter_Junctions
        ;;
    6 | Read_Mapping)
        if [[ ! -d "$GEN_DIR" ]]; then #if genome directory is not specified
            echo "Please specify a valid filepath to the genome directory, exiting..."
            exit 1
        fi
        export RM_PASS="second"
        declare -a files #an array of files
        if [[ -d "$RM_INPUT" ]]; then #if input is a directory
            echo "$RM_INPUT is a directory"
            for f in `find $RM_INPUT -name "*$FORWARD"`; do
                if [[ -f "$f" ]]; then
                    files=("${files[@]}" "$f")
                else
                    echo "$f is not a file"
                fi
            done
            Maxarray=${#files[@]}
        elif [[ -f "$RM_INPUT" ]]; then #if input is a file
            echo "$RM_INPUT is a file"
            Maxarray=$(< $RM_INPUT wc -l)
        else
            echo "Please specify a valid directory or list of input files in the config"
            exit 1
        fi
        if [ "$PE" == "True" ]; then
            echo "$(basename $0): Mapping PE Reads..." >&2
        elif [ "$PE" == "False" ]; then
            echo "$(basename $0): Mapping SE Reads..." >&2
        else
            echo "Please specify in the config file whether data is PE (True/False), exiting..."
            exit 1
        fi
        if [[ -f "$FILTERED_JUNC_LIST" ]]; then #if junction list variable is set as a file
            file_content=$(head -1 ${FILTERED_JUNC_LIST})
            if [[ -f "${file_content}" ]]; then #if first line is a file
                declare -a junctions ### make an array of filtered junction files
                while read line; do
                    junctions=("${junctions[@]}" "$line")
                done < $FILTERED_JUNC_LIST
                export JUNCTIONS="${junctions[@]}"
                export NUM_JUNCTIONS="${#junctions[@]}"
            else #if file input is not a list of files, but list of junctions
                export JUNCTIONS="${FILTERED_JUNC_LIST}"
                export NUM_JUNCTIONS="one"
            fi
            echo "In second-pass mode using ${NUM_JUNCTIONS} junction files"
        elif [[ -z "$FILTERED_JUNC_LIST" ]]; then #if no input here
            echo "Read Mapping without incorporating un-annotated junctions"
        else
            echo "A junction list is specified in config but is not a valid file, exiting..."
            exit 1
        fi
        echo "Max array index is ${Maxarray}">&2
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Read_Mapping.sh" | qsub -l "${RM_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Read_Mapping -V -t 1-"${Maxarray}"
        ;;
    7 | Merge_BAM)
        echo "$(basename $0): Merging BAM files..." >&2
        if [[ -f "$ID_NAMES" ]]; then
            Maxarray=$(cat $ID_NAMES | wc -l)
            echo "Max array index is ${Maxarray}" >&2
            echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Merge_BAM.sh" | qsub -l "${MB_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Merge_BAM -t 1-"${Maxarray}"
        else
            echo "Please specify a valid file containing a list of ID names in the Config file"
        fi
        ;;
    8 | Reference_Prep)
        echo "$(basename $0): Preparing Reference for Quantification..." >&2
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Ref_Prep.sh" | qsub -l "${RP_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Reference_Prep
        ;;
    9 | Transcript_Quant)
        if [ "$PE" == "True" ]; then
            echo "$(basename $0): Quantifying Transcripts for PE data..." >&2
        elif [ "$PE" == "False" ]; then
            echo "$(basename $0): Quantifying Transcripts for SE data..." >&2
        else
            echo "Please specify whether data are paired-end or single-end" >&2
            exit 1
        fi
        for f in $TQ_INPUTDIR/*.bam; do
            if [[ -f "$f" ]]; then
                files=("${files[@]}" "$f")
            else
                echo "$f is not a file"
            fi
        done
        Maxarray=${#files[@]}
        echo "Max array index is ${Maxarray}">&2
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Transcript_Quant.sh" | qsub -l "${TQ_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Transcript_Quant -t 1-"${Maxarray}"
        ;;
    * )
esac