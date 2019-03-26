#PBS -S /bin/bash
#PBS -q batch
#PBS -N Trimm
#PBS -l nodes=1:ppn=4:AMD
#PBS -l walltime=12:00:00
#PBS -l mem=10gb

#PBS -M dittmare@gmail.com
#PBS -m abe

cd $PBS_O_WORKDIR

module load Trimmomatic/0.36-Java-1.8.0_144


INPUTDIR="/scratch/eld72413/Salty_Nut/raw_files/Plate2/Run1" #where the raw sequence data is

OUTPUTDIR="/scratch/eld72413/Salty_Nut/TrimmedReads" #where to put the trimmed sequences

ADAPTERFILE="/home/eld72413/SaltNut/illumina_adapters.txt" #path to adapter file

for d in $INPUTDIR/*ds*; do
if [[ -d "$d" ]]; then
name=$(basename ${d%%-ds.*}"")
for f1 in "$d"/*R1_001.fastq.gz; do
if [[ -f "$f1" ]]; then
f2=${f1%%1_001.fastq.gz}"2_001.fastq.gz"
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -threads 4 -trimlog $OUTPUTDIR/TrimLog $f1 $f2 $OUTPUTDIR/${name}_R1_paired.fq.gz $OUTPUTDIR/${name}_R1_unpaired.fq.gz $OUTPUTDIR/${name}_R2_paired.fq.gz $OUTPUTDIR/${name}_R2_unpaired.fq.gz ILLUMINACLIP:$ADAPTERFILE:2:30:10:1:true LEADING:3 TRAILING:3 MINLEN:20
fi
done
fi
done
