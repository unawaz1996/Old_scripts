#!/bin/bash
# Genetating a script for RNA-seq counts using featurecounts

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -n 12
#SBATCH -N 1
#SBATCH --time=04:00:00
#SBATCH --mem=8GB

# notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=urwah.nawaz@adelaide.edu.au

# loading modules
module load  Subread/1.5.2-foss-2016b

#directories 
#genomeDIR=/data/neurogenetics/RefSeq/GFFandBEDfiles

touch $LOGS/featurecount.txt

INDIR=/fast/users/a1654797/neurogenetics/Debrah/STAR/*.out.bam
OUTDIR=/fast/users/a1654797/neurogenetics/Debrah/Quantification 

if [ -d $OUTDIR ]; then
  echo $OUTDIR "exsits"
else
    mkdir -p $OUTDIR
fi

OUTFILE=$OUTDIR/featurecounts.txt

genomeDIR=/fast/users/a1654797/References
GTF=$genomeDIR/hg19/gencode.v31lift37.annotation.gtf

featureCounts -T 12 -a $GTF -o $OUTFILE $INDIR -s 2 -p
