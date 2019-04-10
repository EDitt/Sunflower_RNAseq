#!/bin/bash

#   Name this project
PROJECT=

#   What email should we use for job notifications?
EMAIL=

############################################
##########   Quality_Assessment   ##########
############################################

#   What are our QSub settings for Quality_Assessment?
#       Below are the recommended settings
QA_QSUB="mem=1gb,nodes=1:ppn=4,walltime=6:00:00"

#   Where are the directories that store your raw data?
#	Sunflower_RNAseq is designed for use on the raw data as obtained from the GGBC
#	There is a directory for each sample/lane run (containing 2 files if paired-end)
#	Directories are named in the following format:
#	SampleID_PlateWell_LaneNum-ds.letters/numbers
#       Include the full file path to these directories containing your raw data
QA_INPUTDIR=""

#   Where do you want the FastQC results?
QA_OUTPUTDIR=""

#	A place for FASTQC to store temporary files
QA_TEMP=""

############################################
##########    Adapter_Trimming    ##########
############################################

#   What are our QSub settings for Adapter_Trimming?
#   Below are the recommended settings
AT_QSUB="mem=10gb,nodes=1:ppn=4,walltime=450:00:00"

#   Where are the directories?
#   Include the full file path to the raw directories
AT_INPUTDIR=""

#   Where do you want the Trimmed Samples to go?
AT_OUTPUTDIR=""

#	What is our adapter file? Include the full file path.
ADAPTERFILE=""

#	Is data paired-end? ("True" or "False")
PE=True

#   What shared suffix do the forward samples have?
#       Example: _1_sequence.txt.gz
FORWARD_NAMING=R1_001.fastq.gz

#   What shared suffix do the reverse samples have?
#       Example: _2_sequence.txt.gz
#	(Only relevant for paired-end data)
REVERSE_NAMING=R2_001.fastq.gz

#	The maximum mismatch count allowed for the "seed" (small section of adapter),
#	which causes the entire alignment between the read and adapter to be scored.
SEEDMISMATCH=2

#	The minimum alignment score threshold for clipping adapter sequence
#	A palindroma approach is used to check for adapter 'read-through'
#	This strategy is only used in PE data, but a value must still be supplied if SE data
PALINDROMECLIP=30

#	The minimum alignment score threshold for clipping adapter sequence
SIMPLECLIP=10

#	The minimum length of adapter sequence to be removed
#	Only relevant for Paired-end data
MINADAPTERLEN=1

#	Whether to keep the reverse read if adapter read-through has been detected by palindrome mode
#	the default behavior is to entirely drop the reverse read
#	Only relevant for Paired-end data
KEEPREADS=true

#	Low quality bases are removed from the beginning of the sequence.
#	What is the minimum quality value required to keep a base at the beginning?
LEADCUT=3

#	Low quality bases are removed from the end of the sequence.
#	What is the minimum quality value required to keep a base at the end?
TRAILCUT=3

#	Reads below a specified minimum length are removed (dropped after other processing steps)
#	What is the minimum length required of reads to be kept?
MINLENGTH=20

############################################
##########      Dependencies      ##########
############################################

#   This section defines installations to
#       various dependencies for Sunflower_RNAseq
#   With module paths specific to the Georgia Advanced Computing Resource Center (GACRC)

#   Do we have FastQC installed?
module load FastQC/0.11.8-Java-1.8.0_144

#	Do we have Trimmomatic installed?
module load Trimmomatic/0.36-Java-1.8.0_144
