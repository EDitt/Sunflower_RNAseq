# Sunflower_RNAseq
A pipeline for analyzing sunflower expression responses to abiotic stress
Inspired by the Morrell lab Sequencing Handling Pipeline: https://github.com/MorrellLAB/sequence_handling

## Programs Used:  
Trimmomatic: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf  
FASTQC: https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf  
STAR: http://chagall.med.cornell.edu/RNASEQcourse/STARmanual.pdf  
RSEM: https://deweylab.github.io/RSEM/README.html

To run `Sunflower_RNAseq`, use the following command, assuming you are in the `Sunflower_RNAseq` directory:
`./Sunflower_RNAseq <handler> Config`
Where `<handler>` is one of the handlers listed below, and `Config` is the full file path to the configuration file

## Pre-processing
Raw data is uploaded/stored in a project folder in the /project/jmblab/ directory
(only accessible through the xfer node)

To analyze this raw data, simply copy this data into your working scratch directory. Do not manipulate the directory structure or names of directories, as this pipeline is designed to be used on the raw data in the format it comes from BaseSpace. 
The raw data comes in directories (of the format "SampleID_PlateWell_LaneNum-ds.letters/numbers"), each containing fastq.gz files (2 files if paired-end sequencing). Manipulating this directory structure or directory names will cause the Sunflower_RNAseq program to not work.


## Step 1: Quality Assessment
Run Quality_Assessment on your raw FastQ files. The Quality_Assessment handler is designed for raw data as it comes from Basespace- i.e. directories for each sample (of the format "SampleID_PlateWell_LaneNum-ds.letters/numbers"), each containing fastq.gz files (2 files if paired-end sequencing).

To run Quality_Assessment, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Quality_Assessment can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)
`./Sunflower_RNAseq Quality_Assessment Config`
where `Config` is the full file path to the configuration file

To get summaries for FastQC results, it is recommended MultiQC is recommended. After Quality_Assessment has finished running, load this module:
`module load MultiQC/1.5-foss-2016b-Python-2.7.14`
And then while in the output directory containing your FASTQC results, simply run
`multiqc .`
This program will then output summary statistics from your FastQ results

## Step 2
Copy data into working 'scratch' directory.

Count the number of files to make sure you have the number you expect
`ls -1 | wc -l`

## Step 3
Use Trimmomatic to trim adapter sequence (see script _**Trimm.sh**_)

## Step 4
Use FASTQC to check quality of data and trimming

## Step 5
Generate genome index for mapping using STAR (only needs to be done once) (see script _**Genome_Index.sh**_)  
You will use the contents of the output file for the next step

## Step 6
Map reads to your genome index using STAR (see script _**Read_Mapping.sh**__)
  - I mapped reads from separate lanes/runs separately - this allows me to test for batch effects after this step and then combine the bam files from the same samples before proceeding

## Step 7
First, prepare the reference for RSEM (see script __**RSEM_prep_ref.sh**__)