#PBS -S /bin/bash
#PBS -q batch
#PBS -N Read_Mapping
#PBS -l nodes=1:ppn=4
#PBS -l walltime=480:00:00
#PBS -l mem=50gb
#PBS -t 0-383%20

#PBS -M brittnein@gmail.com
#PBS -m abe

cd $PBS_O_WORKDIR

module load STAR/2.6.1c-foss-2016b

INPUTDIR="/scratch/bnp34716/Sunflower/Trimmed3" #where the trimmed reads are

OUTPUTDIR="/scratch/bnp34716/Sunflower/Mapped3" #where to put the trimmed sequences

GENOMEDIR="/scratch/bnp34716/Sunflower/GenomeDirNew" #path to Genome Directory

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
