# Sunflower_RNAseq
A pipeline for analyzing sunflower expression responses to abiotic stress

## Step 1
Upload raw data into /project/jmblab/
This folder is backed up and where our highest allotment of storage is

## Step 2
Copy data into working 'scratch' directory.

I did this by creating a list of the filepaths to all folders: `find $(pwd -P) -name "*fastq.gz" | sort -V > sample_list_name.txt`

Then, I copied each file into one new folder (this allows downstream operations to be performed more easily because there will no longer be subdirectory structure to the data).

`mkdir RawData`

`while read line;
do cp $line filepath/to/RawData; done < /filepath/to/sample_list_name.txt`

Count the number of files to make sure you have the number you expect
`ls -1 | wc -l`

## Step 3
Use Trimmomatic to trim adapter sequence

