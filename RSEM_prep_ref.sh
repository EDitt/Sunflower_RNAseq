#PBS -S /bin/bash
#PBS -q batch
#PBS -N Sun_GenomeIndex
#PBS -l nodes=1:ppn=4
#PBS -l walltime=480:00:00
#PBS -l mem=50gb

#PBS -M dittmare@gmail.com
#PBS -m abe

cd $PBS_O_WORKDIR

module load RSEM/1.3.1-foss-2016b
module load STAR/2.6.1c-foss-2016b

#/usr/local/modulefiles_eb/all/RSEM/RSEM/1.3.1-foss-2016b
#--star-path /usr/local/modulefiles_eb/all

rsem-prepare-reference --star  --gtf /scratch/eld72413/XRQ_GenomeDir/XRQ_June2018.gtf /scratch/eld72413/XRQ_GenomeDir/XRQ_June2018.fa /scratch/eld72413/XRQ_GenomeDir/RSEM_ref/XRQ
