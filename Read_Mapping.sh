#PBS -S /bin/bash
#PBS -q batch
#PBS -N Read_Mapping
#PBS -l nodes=1:ppn=4
#PBS -l walltime=480:00:00
#PBS -l mem=50gb
#PBS -t 1-384%20

#PBS -M dittmare@gmail.com
#PBS -m abe

cd $PBS_O_WORKDIR

module load STAR/2.6.1c-foss-2016b

INPUTDIR="/scratch/eld72413/Salty_Nut/TrimmedReads" #where the trimmed reads are

OUTPUTDIR="/scratch/eld72413/Salty_Nut/StarOut4/Plate2" #where to put the trimmed sequences

GENOMEDIR="/scratch/eld72413/XRQ_GenomeDir/GenomeDirNew" #path to Genome Directory

declare -a files
for f in $INPUTDIR/*R1_paired.fq.gz; do
	if [[ -f "$f" ]]; then
		files=("${files[@]}" "$f")
		fi
	done

if [[ -f "${files[${PBS_ARRAYID}]}" ]]; then
	f1="${files[${PBS_ARRAYID}]}"
	f2=${f1%%1_paired.fq.gz}"2_paired.fq.gz"
	name=$(basename ${f1%%_R1_paired.fq.gz}"_")
	echo "Mapping $name Reads"
	/usr/local/apps/eb/STAR/2.6.1c-foss-2016b/bin/STAR \
	--runThreadN 4 \
	--outFileNamePrefix $OUTPUTDIR/"$name" \
	--outFilterScoreMinOverLread 0 \
	--outFilterMatchNminOverLread 0 \
	--outFilterMatchNmin 0 \
	--quantMode TranscriptomeSAM \
	--genomeDir $GENOMEDIR \
	--readFilesCommand gunzip -c \
	--readFilesIn $f1 $f2
fi