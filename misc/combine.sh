#!/bin/bash 

#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -n 12
#SBATCH -N 1
#SBATCH --time=06:00:00
#SBATCH --mem=10GB

# notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=urwah.nawaz@student.adelaide.edu.au


awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' 24-HI1_S19_COMBINED_L000_R1_001.fastq.gz
