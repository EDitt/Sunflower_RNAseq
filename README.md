# Sunflower_RNAseq
A pipeline for analyzing sunflower expression responses to abiotic stress

## Programs Used:  
Trimmomatic: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf  
FASTQC: https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf  
STAR: http://chagall.med.cornell.edu/RNASEQcourse/STARmanual.pdf  
RSEM: https://deweylab.github.io/RSEM/README.html

## Step 1
Upload raw data into /project/jmblab/
This folder is backed up and where our highest allotment of storage is

## Step 2
Copy data into working 'scratch' directory.

I did this by creating a list of the filepaths to all folders: `find $(pwd -P) -name "*fastq.gz" | sort -V > sample_list_name.txt`

Then, I copied each file into one new folder (this allows downstream operations to be performed more easily because there will no longer be subdirectory structure to the data, as there would be if you simply copied the entire folder).

`mkdir RawData`

`while read line;
do cp $line filepath/to/RawData; done < /filepath/to/sample_list_name.txt`

Count the number of files to make sure you have the number you expect
`ls -1 | wc -l`

## Step 3
Use Trimmomatic to trim adapter sequence (see script _**Trimm.sh**_)

Move trimmed, paired reads to a new directory 
`mv /filepath/to/RawData/*_paired.fq.gz /filepath/to/Paired`

Change file extensions to simplify
`for file in *paired.fq.gz; do mv "$file" "${file%_001.fastq.gz_paired.fq.gz}_paired.fq.gz"; done`
This strips the 001.fastq.gz from the filenames

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

#test3

