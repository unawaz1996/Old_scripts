#!/bin/bash

# Got this off Biostars
### unfortunately all fastq files need to be in the same directory for this script to run 

 for i in $(ls -d */ | while read F; do basename $F | rev | cut -c 22- | rev | uniq; done)
    do echo "Merging R1"
    cat "$i"_L00*_R1_001.fastq.gz > "$i"_ME_L000_R1_001.fastq.gz
    echo "Merging R2"
    cat "$i"_L00*_R2_001.fastq.gz > "$i"_ME_L000_R2_001.fastq.gz
done;
