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
	1 | Select_Variants)
        ;;
    2 | Read_Mapping)
        if [[ ! -d "$GEN_DIR" ]]; then #if genome directory is not specified
            echo "Please specify a valid filepath to the genome directory, exiting..."
            exit 1
        fi
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
            echo "$(basename $0): Mapping PE Reads..." >&2
        elif [ "$PE" == "False" ]; then
            echo "$(basename $0): Mapping SE Reads..." >&2
        else
            echo "Please specify in the config file whether data is PE (True/False), exiting..."
            exit 1
        fi
        declare -a junctions ### make an array of filtered junction files
        while read line; do
            junctions=("${junctions[@]}" "$line")
        done < $FILTERED_JUNC_LIST
        export JUNCTIONS="${junctions[@]}"
        export NUM_JUNCTIONS="${#junctions[@]}"
        echo "Variant-Aware read mapping using ${NUM_JUNCTIONS} junction files"
        echo "Max array index is ${Maxarray}">&2
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/VariantAware_Read_Mapping.sh" | qsub -l "${RM_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Read_Mapping  -V -t 1-"${Maxarray}"
        ;;
    3 | Read_Counter)
		;;
	* )
esac