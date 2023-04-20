#!/bin/bash

#### Test run for aligning reads to graph genome and calling variants 
#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -n 32
#SBATCH -N 1
#SBATCH --time=05:00:00
#SBATCH --mem=80GB

# notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=urwah.nawaz@adelaide.edu.au#!/bin/bash

## Loading modules 
module load vg/1.10.0 


### directions and files 

genome=/data/neurogenetics/graphPangenome/pangenome
#INDIR=/uofaresstor/neurogenetics/sequences/Illumina/genome/AGRF_KCCG_WGS
#INDIR=/fast/users/a1654797/neurogenetics/graphGenome/sequences
INDIR=/fast/users/a1654797/neurogenetics/RNA_seq_analysis/DATA/Deepti_mRNA-Seq_Aug2018-91474384/Combied_runs


vg map -f ${INDIR}/Control01_R1.fastq -f ${INDIR}/Control01_R2.fastq -x ${genome}/1kGP_SGDP_Hominins.xg -g ${genome}/1kGP_SGDP_Hominins.gcsa -t 32 -s 0 -u 1 -m 1 > aln.gam
