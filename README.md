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
`./Sunflower_RNAseq.sh Genome_Index Config`  
where `Config` is the full file path to the configuration file.

You will use the contents of the output (directory specified in the Config file) for the next step

## Step 4: Read_Mapping

The Read_Mapping handler uses STAR to map reads to the genome indexed in step 3.

This handler can accept as input EITHER a directory or a text-file list of forward samples (it will find the reverse samples based on the naming suffix specified in the config file)

### Option for 2-pass mapping
STAR can perform a 2-pass mapping strategy to increase mapping sensitivity around novel splice junctions. This works by running a 1st mapping pass for all samples with the "usual" parameters to identify un-annotated splice junctions. Then mapping is re-run, using the junctions detected in the first pass (in addition to the annotated junctions specified in the genome annotation file). This 2-pass mapping strategy is recommended by GATK and ENCODE best-practices for better alignments around novel splice junctions.

While STAR can perform 2-pass mapping on a per-sample basis, in a study with multiple samples, it is recommended to collect 1st pass junctions from all samples for the second pass. Therefore, the recommended 2-pass procedure is described below: 

#### Step 4a: Collect_Junctions
This step identifies novel junctions during a first read-mapping pass and outputs them as "SJ.out.tab" files for each sample. In the handler used here, only junctions supported by at least 1 uniquely mapped read will be output. This step shares variables for Read_Mapping in the config file, so make sure these are filled out.

Once the variables have been defined, Collect_Junctions can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)  
`./Sunflower_RNAseq.sh Collect_Junctions Config`   
where `Config` is the full file path to the configuration file.

#### Step 4b: Filter_Junctions
This step will concatenate the junction files discovered in Step 4a across samples, and then filter them based on user-defined parameters in the config file. This step is not required, as STAR will automatically concatenate the "SJ.out.tab" files before mapping. If you want to skip this step, you can instead just pass a list of all the "SJ.out.tab" files for all samples from Step 4a directly into the `FILTERED_JUNC_LIST` variable for Read_Mapping. (In fact, this is the 2-pass mapping procedure describe in the STAR manual). However, filtering junctions is recommended when you have large numbers of samples, and there are several reasons we have implemented this intermediate filtering step: 
1.) Spurious junctions may increase the number of reads that map to multiple places in the genome   
2.) Filtering and concatenating junction files before mapping speeds up the mapping step in 4c - both because of the smaller number of junctions and because STAR then doesn't need to perform the concatenation across large numbers of samples for each sample separately.  
3.) These lists can be more easily saved in case one wants to redo the mapping  

Once the variables in the configuration file have been defined, Filter_Junctions can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)  
`./Sunflower_RNAseq.sh Filter_Junctions Config`   
where `Config` is the full file path to the configuration file.

#### Step 4c: Read_Mapping  
This step will take all of the novel junctions discovered in the first pass and use them to re-map reads for each sample (in addition to the already-annotated junctions from your annotation file). A .txt file with a list of all junction files needs to be specified in the config file for the `FILTERED_JUNC_LIST` variable. As mentioned previously, this can be a list of un-filtered, un-concatenated "SJ.out.tab" files from all samples (if skipping step 4b) or the names of the filtered junction lists output from step 4b. 

Because we often process samples in batches (based on flow cell run), we have multiple outputs from step 4b, and simply make a list of these as input for the `FILTERED_JUNC_LIST` variable. As mentioned previously, STAR can handle multiple junction files as input (and will concatenate before mapping).

To run Read_Mapping, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Read_Mapping can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)  
`./Sunflower_RNAseq.sh Read_Mapping Config`   
where `Config` is the full file path to the configuration file.

### Option for 1-pass mapping  
If you want to map your reads without the addition of novel junctions discovered in a first mapping step, you can skip the "Collect_Junctions" and "Filter_Junctions" steps and leave the `FILTERED_JUNC_LIST` variable blank. All other variables need to be specified for Read_Mapping in the configuration file. Once the variables have been defined, Read_Mapping can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)  
`./Sunflower_RNAseq.sh Read_Mapping Config`   
where `Config` is the full file path to the configuration file.

### Final Notes:  
#### Read groups  
If you have sequence data from the same sample across multiple lanes/runs, the best practice is to map these separately (in order to test for batch effects), and then combine resulting bam files (step 5) for each sample before proceeding to transcript quantification. The read mapping step will add read groups (@RG) headers during mapping based upon file names and variables specified in the config file to differentiate separate read groups after merging.

#### Alignment File Output  
Two alignment files are output from STAR with this handler - one alignment file in genomic coordinates and one translated into transcript coordinates (e.g.`*Aligned.toTranscriptome.out.bam`). The latter is needed for transcript quantification with RSEM (step 7). The default format of the former is as an unsorted SAM file. If you plan on using the genomic coordinate alignments for SNP calling, you have the option of getting these output as coordinate-sorted BAM files (similar to the `samtools sort` command) by putting a "yes" for the `GENOMIC_COORDINATE_BAMSORTED` variable in the config. Note that this will add significant computational time and memory.

## Step 5: Merge_BAM (Optional)

If you choose to map reads from different runs/lanes for the same sample (recommended), the Merge_BAM handler will merge the BAM files from the STAR output (the "Aligned.toTranscriptome.out.bam" files) using samtools before proceeding to transcript quantification.

In addition to specifying input and output directories where these files are located, this handler requires a .txt file listing all of the sample ID names. These should correspond to the leading sample ID name in your input BAM files.

Example: the sample ID for files 145-2068-L001_Aligned.toTranscriptome.out.bam and 145-2068-2_L003_Aligned.toTranscriptome.out.bam is 145.

To run Merge_BAM, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Merge_BAM can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)
`./Sunflower_RNAseq.sh Merge_BAM Config`

## Step 6: Reference_Prep

The Reference_Prep handler uses RSEM to prepare reference transcripts used for transcript quantification

RSEM can extract reference sequences from a genome if it is provided with gene annotations in a GTF/GFF3 file. If the annotation file is in GFF3 format, RSEM will first convert it to GTF format with the file name 'reference_name.gtf'

Alternatively, you can provide RSEM with transcript sequences directly in the form of fasta files. To do this the `--gtf` or `--gff3` flags should be commented out of the Ref_Prep.sh script. In this case, RSEM assumes the reference fasta files contain the reference transcripts and that the name of each sequence in the Multi-FASTA files are transcript IDs.

To run Reference_Prep, all common and handler-specific variables must be defined within the configuration file. Once the variables have been defined, Reference_Prep can be submitted to a job scheduler with the following command (assuming that you are in the directory containing `Sunflower_RNAseq`)
`./Sunflower_RNAseq.sh Reference_Prep Config`

## Step 7: Transcript_Quant

The Transcript_Quant handler uses RSEM to calculate expression using the reference prepared in the previous "Reference Prep" step. When running this handler, make sure that the variables `RSEM_ref` and `REF_NAME` are supplied (under Reference_Prep) in addition to the handler-specific variables. Sunflower_RNAseq uses RSEM's `--bam` option to allow input files to be in BAM format. The Transcript_Quant handler assumes that data are from STAR alignments.

If data are from single-end reads, providing a fragment-length mean and fragment-length standard deviation is important for the accuracy of expression levels and you should fill out these variables in the config file. If data are paired-end, these variables are ignored (RSEM will automatically learn a fragment length distribution from the data).

## Step 8: Make a data matrix

Once you have output from the transcript quantification step, you can make a data matrix of your expression results using the `rsem-generate-data-matrix` script. Example code is below: 

`module load RSEM/1.3.1-foss-2016b`  
`rsem-generate-data-matrix InputDir/*.genes.results > OutputDir/name.matrix` 

You can now load this matrix into R to begin differential expression analysis!
