library(tibble)
library(pander)
library(kableExtra)
library(dplyr)
library(magrittr)


#### loading the data files 

gtex_files = read.table("files//gtex_runNames.txt", header=FALSE)
SRA_info= read.table("files/SraRunTable.txt", header=TRUE, sep = "\t")

#### comparing the files 
gtex_files$V1 <- gsub(".sra", "", gtex_files$V1)
gtex_files
dim(SRA_info)
head(SRA_info)


matched_files <- SRA_info[SRA_info$Run %in% gtex_files$V1,]
head(matched_files)


rna_seq <- matched_files %>% 
  dplyr::filter(molecular_data_type == "RNA Seq (NGS)")
rna_seq

dim(rna_seq)


length(unique(rna_seq$submitted_subject_id))
