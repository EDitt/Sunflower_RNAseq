#!/bin/bash

#   Name this project
PROJECT=

#   What email should we use for job notifications?
EMAIL=

#	Is data paired-end? ("True" or "False")
PE=True


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