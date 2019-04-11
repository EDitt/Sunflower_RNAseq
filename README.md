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
The raw data comes in directories (of the format "SampleID_PlateWell_LaneNum-ds.letters/numbers"), each containing fastq.gz files (2 files if paired-end sequencing). Manipulating this directory structure or directory names will cause the Sunflower_RNAseq program to not work properly.


## Step 1: Quality Assessment
Run Quality_Assessment on your raw FastQ files. The Quality_Assessment handler is designed for raw data as it comes from Basespace- i.e. directories for each sample (of the format "SampleID_PlateWell_LaneNum-ds.letters/numbers"), each containing fastq.gz files (2 files if paired-end sequencing).

To run Quality_Assessment, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Quality_Assessment can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)
`./Sunflower_RNAseq Quality_Assessment Config`
where `Config` is the full file path to the configuration file

To get summaries for FastQC results, MultiQC is recommended. After Quality_Assessment has finished running, load this module:
`module load MultiQC/1.5-foss-2016b-Python-2.7.14`
And then while in the output directory containing your FASTQC results, simply run
`multiqc .`
This program will then output summary statistics from your FastQ results

## Step 2: Adapter Trimming
The Adapter_Trimming handler uses Trimmomatic to trim adapter sequences from FastQ files. Trimmomatic takes paired-end information into account when doing so (if applicable). This handler uses the raw fastq.gz files downloaded directly from Basespace (Directories named: "SampleID_PlateWell_LaneNum-ds.letters/numbers"). This handler will work with paired-end or single-end sequencing data.

To run Adapter_Trimming, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Adapter_Trimming can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)
`./Sunflower_RNAseq Adapter_Trimming Config`
where `Config` is the full file path to the configuration file

In addition to adapter trimming, Trimmomatic can also perform quality trimming. However, the Adapter_Trimming handler used here does not use Trimmomatic's quality trimming options. Many caution against quality trimming, as it is believed to be unnecessary since read mapping approaches can take quality scores into account. If you do want to use Trimmomatic's quality trimming capabilities, the `Trimm.sh` code must be modified and new variables defined in the configuration file. Read the Trimmomatic manual for more information.

It is recommended that you re-run Quality_Assessment after adapter trimming to ensure that any adapter contamination was eliminated.

## Step 3: Generate a Genome Index

This handler will generate a genome index using FASTA and GFF3 formatted annotations. This step only needs to be performed only once for each genome/annotation combination.

If using a GTF file instead of a GFF3 file for genome indexing (the default file format for this step), the line `--sjdbGTFtagExonParentTranscript Parent` should be commented out of the Genome_Index.sh script.

To run Genome_Index, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Genome_Index can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)
`./Sunflower_RNAseq Genome_Index Config` where `Config` is the full file path to the configuration file.

You will use the contents of the output (directory specified in the Config file) for the next step

## Step 4: Read Mapping
Map reads to your genome index using STAR (see script _**Read_Mapping.sh**__)
  - I mapped reads from separate lanes/runs separately - this allows me to test for batch effects after this step and then combine the bam files from the same samples before proceeding

## Step 5: Counting Reads
First, prepare the reference for RSEM (see script __**RSEM_prep_ref.sh**__)