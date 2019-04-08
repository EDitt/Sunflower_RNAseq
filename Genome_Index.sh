#PBS -S /bin/bash
#PBS -q batch
#PBS -N GenomeIndex
#PBS -l nodes=1:ppn=4
#PBS -l walltime=480:00:00
#PBS -l mem=50gb

#PBS -M dittmare@gmail.com
#PBS -m abe

cd /scratch/eld72413/Salty_Nut

module load STAR/2.6.1c-foss-2016b

/usr/local/apps/eb/STAR/2.6.1c-foss-2016b/bin/STAR \
--runThreadN 4 \
--runMode genomeGenerate \
--genomeDir /scratch/eld72413/HA412_GenomeDir \
--genomeFastaFiles /scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.fasta \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFfile /scratch/eld72413/Ha412HOv2.0/Ha412HOv2.0-20181130.gff3 \
--sjdbOverhang 74
