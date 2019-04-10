#!/bin/bash

#   Name this project
PROJECT=Test

#   What email should we use for job notifications?
EMAIL=dittmare@gmail.com

#   What reference genome are we using?
#       Include the full file path.
REF_GEN=

#	Where is the adapter file?
#		Include the full file path
ADAPTERFILE="/home/eld72413/SaltNut/illumina_adapters.txt"

############################################
##########   Quality_Assessment   ##########
############################################

#   What are our QSub settings for Quality_Assessment?
#       Below are the recommended settings
QA_QSUB="mem=1gb,nodes=1:ppn=4,walltime=6:00:00"

#   Where are the directories?
#       Include the full file path to the raw directories
QA_INPUTDIR="/home/eld72413/SaltNut/ScriptTest/Input"

#   Where are we storing the output files?
QA_OUTPUTDIR="/home/eld72413/SaltNut/ScriptTest/Output"

#	A place for FASTQC to store temporary files
QA_TEMP="/scratch/eld72413/Tmp"

############################################
##########      Dependencies      ##########
############################################

#   This section defines installations to
#       various dependencies for Sunflower_RNAseq
#   With module paths specific to the Georgia Advanced Computing Resource Center (GACRC)

#   Do we have FastQC installed?
module load FastQC/0.11.8-Java-1.8.0_144
