#PBS -S /bin/bash
#PBS -q batch
#PBS -N Star_Jan8
#PBS -l nodes=1:ppn=4
#PBS -l walltime=480:00:00
#PBS -l mem=50gb
#PBS -t 1-768%50

#PBS -M dittmare@gmail.com
#PBS -m abe

cd $PBS_O_WORKDIR

module load STAR/2.6.1c-foss-2016b

for f1 in /scratch/eld72413/Salty_Nut/TrimmomaticReads/Paired/*R1_paired.fq.gz
do
	f2=${f1%%1_paired.fq.gz}"2_paired.fq.gz"
	name=$(basename ${f1%%_R1_paired.fq.gz}"_")
	/usr/local/apps/eb/STAR/2.6.1c-foss-2016b/bin/STAR --runThreadN 4 --outFileNamePrefix /scratch/eld72413/Salty_Nut/StarOut4/"$name" --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 --quantMode TranscriptomeSAM --genomeDir /scratch/eld72413/XRQ_GenomeDir/GenomeDirNew --readFilesCommand gunzip -c --readFilesIn $f1 $f2
done
