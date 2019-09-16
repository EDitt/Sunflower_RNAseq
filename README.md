# Sunflower_RNAseq
A pipeline for analyzing sunflower expression responses to abiotic stress  
Inspired by the Morrell lab Sequencing Handling Pipeline: https://github.com/MorrellLAB/sequence_handling

## Programs Used:  
FASTQC: https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf  
MultiQC: https://multiqc.info/  
Trimmomatic: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf  
STAR: http://chagall.med.cornell.edu/RNASEQcourse/STARmanual.pdf  
RSEM: https://deweylab.github.io/RSEM/README.html

To run `Sunflower_RNAseq`, use the following command, assuming you are in the `Sunflower_RNAseq` directory:  
`./Sunflower_RNAseq.sh <handler> Config`  
Where `<handler>` is one of the handlers listed below, and `Config` is the full file path to the configuration file

## Pre-processing
Raw data is uploaded/stored in a project folder in the /project/jmblab/ directory
(only accessible through the xfer node)

To analyze this raw data, simply copy this data into your working scratch directory.

## Step 1: Quality_Assessment
Run Quality_Assessment on your raw FastQ files. 

To run Quality_Assessment, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Quality_Assessment can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)
`./Sunflower_RNAseq.sh Quality_Assessment Config`
where `Config` is the full file path to the configuration file

A directory containing your files to be analyzed must be specified in the config file. It is ok if this is a directory containing sub-directories for each sample (which is the format for raw data as it comes from Basespace).

After quality has been assessed for each sample, the FastQC results will be summarized using MultiQC. These summary statistics will be located in the output directory specified in the config file.

## Step 2: Adapter_Trimming
The Adapter_Trimming handler uses Trimmomatic to trim adapter sequences from FastQ files. Trimmomatic takes paired-end information into account when doing so (if applicable).

The Adapter_Trimming handler can accept as input EITHER a directory (which can be to multiple sub-directories for each sample) or a text-file list of forward samples (it will find the reverse samples based on the naming suffix specified in the config file)

To run Adapter_Trimming, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Adapter_Trimming can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)  
`./Sunflower_RNAseq.sh Adapter_Trimming Config`  
where `Config` is the full file path to the configuration file

While Trimmomatic can also perform quality trimming, the Adapter_Trimming handler used here does not use Trimmomatic's quality trimming options. Many caution against quality trimming, as it is believed to be unnecessary since read mapping approaches can take quality scores into account. If you do want to use Trimmomatic's quality trimming capabilities, the `Trimm.sh` code must be modified and new variables defined in the configuration file. Read the Trimmomatic manual for more information: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf  

It is recommended that you re-run Quality_Assessment after adapter trimming to ensure that any adapter contamination was eliminated.

## Step 3: Genome_Index  

This handler will generate a genome index using FASTA and GFF3 or GTF formatted annotations. This step only needs to be performed once for each genome/annotation combination.

If using a GFF3 file for genome indexing rather than the default GTF file, the option `--sjdbGTFtagExonParentTranscript Parent` is added to the script 

To run Genome_Index, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Genome_Index can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)
`./Sunflower_RNAseq.sh Genome_Index Config` where `Config` is the full file path to the configuration file.

You will use the contents of the output (directory specified in the Config file) for the next step

## Step 4: Read_Mapping

The Read_Mapping handler uses STAR to map reads to the genome indexed in step 3.

To run Read_Mapping, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Read_Mapping can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)
`./Sunflower_RNAseq.sh Read_Mapping Config` where `Config` is the full file path to the configuration file.

If you have sequence data from the same sample across multiple lanes/runs, the best practice is to map these separately (in order to test for batch effects), and then combine resulting bam files for each sample before proceeding to transcript quantification.

## Step 5: Merge_BAM (Optional)

If you choose to map reads from different runs/lanes for the same sample (recommended), the Merge_BAM handler will merge the BAM files from the STAR output before proceeding to transcript quantification.

## Step 6: Ref_Prep

The Ref_Prep handler uses RSEM to prepare reference transcripts used for transcript quantification

RSEM can extract reference sequences from a genome if it is provided with gene annotations in a GTF/GFF3 file. If the annotation file is in GFF3 format, RSEM will first convert it to GTF format with the file name 'reference_name.gtf'

Alternatively, you can provide RSEM with transcript sequences directly in the form of fasta files. To do this the `--gtf` or `--gff3` flags should be commented out of the Ref_Prep.sh script. In this case, RSEM assumes the reference fasta files contain the reference transcripts and that the name of each sequence in the Multi-FASTA files are transcript IDs.

## Step 7: Transcript_Quant

The Transcript_Quant handler uses RSEM to calculate expression