#PBS -S /bin/bash
#PBS -q batch
#PBS -N Trimm_Jan2
#PBS -l nodes=1:ppn=4:AMD
#PBS -l walltime=480:00:00
#PBS -l mem=10gb

#PBS -M dittmare@gmail.com
#PBS -m abe

cd $PBS_O_WORKDIR

module load Trimmomatic/0.36-Java-1.8.0_144


for f1 in /scratch/eld72413/Salty_Nut/raw_files/all/*R1_001.fastq.gz
do 
	f2=${f1%%1_001.fastq.gz}"2_001.fastq.gz"
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 4 $f1 $f2 ${f1}_paired.fq.gz ${f1}_unpaired.fq.gz ${f2}_paired.fq.gz ${f2}_unpaired.fq.gz ILLUMINACLIP:/home/eld72413/SaltNut/illumina_adapters.txt:2:30:10 LEADING:10 TRAILING:10 MINLEN:40
done
