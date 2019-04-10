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
#       Below are the recommended settings
AT_QSUB="mem=1gb,nodes=1:ppn=4,walltime=50:00:00"

#   Where are the directories?
#       Include the full file path to the raw directories
AT_INPUTDIR=""

#   Where do you want the Trimmed Samples to go?
AT_OUTPUTDIR=""

#	What is our adapter file? Include the full file path.
ADAPTERFILE=""

#	Where 
LEADCUT=

#	Where 
TRAILCUT=

#	Where 
MINLEN=

#   What shared suffix do the forward samples have?
#       Example: _1_sequence.txt.gz
FORWARD_NAMING=_R1.fastq.gz

#   What shared suffix do the reverse samples have?
#       Example: _2_sequence.txt.gz
REVERSE_NAMING=_R2.fastq.gz

#   TO DO: Handle single-end reads

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
