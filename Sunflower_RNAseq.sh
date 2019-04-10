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


#   Run sequence_handling
case "${ROUTINE}" in
    1 | Quality_Assessment)
        echo "$(basename $0): Assessing quality..." >&2
        declare -a files #an array of directories to each sample
        for d in $QA_INPUTDIR/*ds*; do
            if [[ -d "$d" ]]; then
                files=("${files[@]}" "$d")
            else 
                echo "Please specify a path to valid directories in the config file"
            fi
        done
        Maxarray=${#files[@]}
        echo "Max array index is ${Maxarray}">&2
        #source "${SUNFLOWER_RNA}"/FASTQC.sh
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/FASTQC.sh" | qsub -l "${QA_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Quality_Assessment -t 1-"${Maxarray}"%20
        ;;
    2 | Adapter_Trimming)
        echo "$(basename $0): Trimming Adapters..." >&2
        declare -a files #an array of directories to each sample
        for d in $AT_INPUTDIR/*ds*; do
            if [[ -d "$d" ]]; then
                files=("${files[@]}" "$d")
            else
                echo "Please specify a path to valid directories in the config file"
            fi
        done
        Maxarray=${#files[@]}
        echo "Max array index is ${Maxarray}">&2
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Trimm.sh" | qsub -l "${AT_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Adapter_Trimming -t 1-"${Maxarray}"%20
        ;;
* )
esac