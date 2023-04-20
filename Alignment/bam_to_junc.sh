#!/bin/bash


#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=01:00:00
#SBATCH --mem=12GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=urwah.nawaz@adelaide.edu.au


### loading modules
module load  SAMtools/1.8-foss-2016b
module load Python/3.6.1-foss-2016b



## scripts, directories and files
leafcutter=/fast/users/a1654797/leafcutter
INDIR=/fast/users/a1654797/neurogenetics/RNA_seq_analysis/Genome_align/hg19
OUTDIR=/fast/users/a1654797/neurogenetics/RNA_seq_analysis/Genome_align/hg19/leafcutter

if [ -d $OUTDIR]; then
  echo $OUTDIR "exists"
else
    mkdir -p $OUTDIR
fi

cd $INDIR

BAMS=$(ls *.bam)

for bamfile in ${BAMS} ;
do
  bam=${file%.bam}
    echo "Converting $bamfile to $bamfile.junc"
    bash $leafcutter/scripts/bam2junc.sh ${INDIR}/${bam}.bam ${OUTDIR}/${bam}.junc
    echo ${bam}.junc >> test_juncfiles.txt
done
