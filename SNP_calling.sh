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
    1 | Merge_BAM)
        echo "$(basename $0): Merging BAM files..." >&2
        if [[ -f "$ID_NAMES" ]]; then
            Maxarray=$(cat $ID_NAMES | wc -l)
            echo "Max array index is ${Maxarray}" >&2
            echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Merge_BAM_Picard.sh" | qsub -l "${MB_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Merge_BAM -t 1-"${Maxarray}"
        else
            echo "Please specify a valid file containing a list of ID names in the Config file"
        fi
        ;;
    2 | Process_BAM)
        echo "$(basename $0): Processing BAM files..." >&2
        declare -a files #an array of files
        for f in `find $PB_INPUTDIR -name "*$PB_SUFFIX"`; do
            if [[ -f "$f" ]]; then
                files=("${files[@]}" "$f")
            else
                echo "$f is not a file"
            fi
        done
        Maxarray=${#files[@]}
        echo "Max array index is ${Maxarray}">&2
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Process_BAM.sh" | qsub -l "${PB_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Process_BAM -t 1-"${Maxarray}"
        ;;
    3 | Validate_BAM)
        echo "$(basename $0): Validating BAM files..." >&2
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Validate_BAM.sh" | qsub -l "${VB_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Validate_BAM
        ;;
	* )
esac