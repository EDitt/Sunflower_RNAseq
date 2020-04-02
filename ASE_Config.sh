#!/bin/bash

#   Name this project
PROJECT=

#   What email should we use for job notifications?
EMAIL=

#	Is data paired-end? ("True" or "False")
PE=True

############################################
##########    Read_Mapping        ##########
############################################

#####IMPORTANT:
# Make sure that the "GEN_DIR" variable is supplied above (under Genome Index)
# This directory will also be used for read mapping

#	What are our QSub settings for Read_Mapping?
#	Below are the recommended settings
#	Note that STAR needs at least 30gb of memory to run
RM_QSUB="mem=50gb,nodes=1:ppn=16,walltime=450:00:00"

#   Where are the trimmed files? This handler will accept EITHER a text file list of forward samples to map (include the full file path for each sample) OR
#	a directory of input files (include the full file path). This may be the same as your output directory defined in adapter trimming (AT_OUTPUTDIR)
RM_INPUT=""

#   What shared suffix do the forward samples have?
#       Example: _R1_paired.fq.gz
FORWARD="_R1_paired.fq.gz"

#   What shared suffix do the reverse samples have? (leave blank if single-end)
#       Example: _R2_paired.fq.gz
REVERSE="_R2_paired.fq.gz"

#   Where do you want your output files to go? Include the full file path
RM_OUTPUTDIR=""

#	Define the number of threads to be used. 
#	This must be set to the number of available cores on the server node
RM_NTHREAD=16

#	What is the maximum number of mismatches allowed?
#	Alignment will output only if it has no more mismatches than this value
MAX_MIS=10

#	What is the maximum number of loci the read is allowed to map to?
#	All alignments will output only if the read maps to no more than this value
#	If greater than this value, read counted as "mapped to too many loci"
MAX_N=10

#	Do you want unmapped or partially mapped reads to output in separate fasta/fastq files?
#	If not, put "None", if yes, "Fastx"
UNMAP_F=Fastx

#	What is the minimum score needed for alignment?
#	Normalized by read length- (sum of mates' lengths for paired-end reads)
MINSCORE_READL=0.66

#	What is the minimum number of matched bases needed for alignment?
#	Normalized by read length- (sum of mates' lengths for paired-end reads)
MINMATCH_READL=0.66

#####	For adding read group (@RG) tags to mapped sequence:
#	This information is useful for differentiating reads coming from different runs/lanes even after merging

#	Name of Flowcell (or other name that differentiates a particular run, e.g. "Run1")
FLOWCELL_NAME=""

# sequencing platform used
PLATFORM="ILLUMINA"


### Variables not shared with Collect_Junctions handler:

#	STAR outputs alignment in genome coordinates as well as transcript coordinates. The latter is used for
#	RSEM (step 7). If you plan on using the genomic coordinate alignment files for SNP calling, you will need 
#	sorted BAM files. If you want STAR to output genome alignments as sorted BAM files, put "yes" here. Note that this
#	increases the time it takes for read mapping to run. If "no", the genome alignment file will be an unsorted SAM.
GENOMIC_COORDINATE_BAMSORTED="no"

#	The type of quantification requested
#	Either "TranscriptomeSAM" -outputs SAM/BAM alignments to transcriptome in a separate file
#	Or "GeneCounts" - read counts per gene
#	Can also put "-" for none
QUANT=TranscriptomeSAM

#	2nd-pass read mapping will use junctions discovered in the first pass (Collect_Junctions) 
#	Need as input a list of all filtered junction files from samples. Include full filepath
#	If not using additional junctions, leave blank
FILTERED_JUNC_LIST=


############################################
######   Variant-Aware Read_Mapping    #####
############################################

#   Genome Directory
#	Include the full filepath
GEN_DIR=""

#####	VARIANT-AWARE MAPPING
# To perform variant-aware approaches, you need to supply a VCF file as input
# This will automatically perform WASP filtering on reads
VARIANT-AWARE="yes"

# Full filepath to VCF
VCF=""