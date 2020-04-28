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
    4 | Haplotype_Caller)
        echo "$(basename $0): Identifying variants using GATK's Haplotype Caller..." >&2
        if [[ -f "$HC_INPUT" ]]; then
            Maxarray=$(cat $HC_INPUT | wc -l)
            echo "Max array index is ${Maxarray}" >&2
            echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Haplotype_Caller.sh" | qsub  -q "${HC_QUEUE}" -l "${HC_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Haplotype_Caller -t 1-"${Maxarray}"
        else
            echo "Please specify a valid input file"
        fi
        ;;
    4 | Genomics_DB_Import)
        echo "$(basename $0): Running Genomics DB Import..." >&2
        if [[ -d "$GD_INPUT" ]]; then #if input is a directory
            echo "$GD_INPUT is a directory"
            declare -a samples #an array of files
            for f1 in `find $GD_INPUT -name "*g.vcf"`; do
                if [[ -f "$f1" ]]; then
                    samples=("${samples[@]}" "-V $f1")
                else
                    echo "Please specify a path to valid files in the config file"
                fi
            done
        elif [[ -f "$GD_INPUT" ]]; then #if input is a file
            echo "$GD_INPUT is a file"
            declare -a samples
            while read line; do
                samples=("${samples[@]}" "-V $line")
            done < $GD_INPUT
        else
            echo "Please specify a valid directory or list of input files in the config"
        fi
        export INPUT="${samples[*]}"
        int_ext=$(basename ${GD_INTERVALS##*.}) # is interval file .bed format?
        if [[ "$int_ext" == "bed" ]]; then # if interval extension is .bed convert to Picard interval list
            GEN_DICT="${GEN_FASTA%%fasta}dict"
            if [[ -f "$GEN_DICT" ]]; then # if reference dictionary exists
                echo "Converting .bed file interval format to 1-based coordinate system"
                int_name=$(basename ${GD_INTERVALS%%.bed})
                java -jar ${PICARD_JAR} BedToIntervalList \
                I=${GD_INTERVALS} \
                O=${GD_OUTPUTDIR}/${int_name}.interval_list \
                SD=${GEN_DICT}
                export IntList="${GD_OUTPUTDIR}/${int_name}.interval_list"
            else
                echo "No sequence dictionary file (.dict) detected in the $(dirname "${GEN_FASTA}") directory"
                echo "Please generate this file and re-run, exiting..."
                exit 1
            fi
        else
            echo "Detected 1-based coordinate interval list"
            export IntList="${GD_INTERVALS}"
            #IntNum=$(< $IntList wc -l)
        fi
        IntNum=$(grep -v "^@" $IntList | wc -l)
        if [[ "$GD_SCAFFOLDS" != "false" ]]; then #if scaffold sequence added
            echo "Adding scaffold regions to process in addition to $IntNum chromosomal regions"
            export Maxarray=$(($IntNum + 1))
        else
            echo "No scaffold regions specified in config"
            Maxarray="$IntNum"
        fi
        echo "Max array index is ${Maxarray}" >&2
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Genomics_DB_Import.sh" | qsub  -q "${GD_QUEUE}" -l "${GD_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Genomics_DB_Import -V -t 1-"${Maxarray}"
        ;;
    5 | Genotype_GVCFs)
        echo "$(basename $0): Genotyping GVCFs using workspace created by Genomics_DB_Import..." >&2
        int_ext=$(basename ${GD_INTERVALS##*.})
        if [[ "$int_ext" == "bed" ]]; then # if interval extension is .bed then use created Picard interval list
            echo "Using Picard interval list created in previous step"
            int_name=$(basename ${GD_INTERVALS%%.bed})
            export IntList="${GD_OUTPUTDIR}/${int_name}.interval_list"
        else
            echo "Using input interval list"
            export IntList="${GD_INTERVALS}"
        fi
        IntNum=$(grep -v "^@" $IntList | wc -l)
        if [[ "$GD_SCAFFOLDS" != "false" ]]; then #if scaffold sequence added
            echo "Adding scaffold regions to process in addition to $IntNum chromosomal regions"
            export Maxarray=$(($IntNum + 1))
        else
            echo "No scaffold regions specified in config"
            Maxarray="$IntNum"
        fi
        echo "Max array index is ${Maxarray}" >&2
        echo "source ${CONFIG} && source ${SUNFLOWER_RNASEQ}/Genotype_GVCFs.sh" | qsub  -l "${GV_QSUB}" -e "${ERROR}" -o "${ERROR}" -m abe -M "${EMAIL}" -N "${PROJECT}"_Genotype_GVCFs -V -t 1-"${Maxarray}"
        ;;
    6 | Combine_VCFs)
        echo "$(basename $0): Combining VCFs and filtering variants..." >&2
        ;;
	* )
esac