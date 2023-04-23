### creating the pdata files and making sure 
### everything is labelled appropriately 


library(ngsReports)
library(dplyr)
library(magrittr)
library(pander)
library(kableExtra)
library(tibble)
library(knitr)
library(ggplot2)
library(summarytools)


sraData <- read.table("files/SraRunTable.txt", sep = "\t", header=TRUE)  ## the data here is SRA data attributes- Run is unique per line 
phenotype<-read.table("files/phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt",sep="\t",
                      skip=10,header=T,fill=TRUE,quote="") ### phenotype data - SUBJID is unique 



fileDir <- ("raw_fastqc")
files <- list.files(fileDir, pattern = "fastqc.zip$", full.names = TRUE)
fdl <- getFastqcData(files)
total_files <- gsub("raw_fastqc/","", files)
total_files <- gsub("_1_fastqc.zip","", total_files)
total_files <- gsub("_2_fastqc.zip","", total_files)


sra_inData <- sraData[sraData$Run %in% total_files,]
sra_inData_rnaseq <- sra_inData %>% 
  dplyr::filter(molecular_data_type == "RNA Seq (NGS)", Center_Name == "BI")  %>% 
  as.data.frame()
sra_inData_rnaseq


### can confirm is unique - not sure why the subject ID does not work
subject_to_sample_inData <- subject_to_sample[subject_to_sample$BioSample.Accession %in% sra_inData_rnaseq$BioSample,]


### phenotype data - can confirm this is correct  
phenotypeData <- phenotype[phenotype$dbGaP_Subject_ID %in% subject_to_sample_inData$dbGaP_Subject_ID,]


sample_attributes_inData <- sample_attributes[sample_attributes$dbGaP_Sample_ID %in% subject_to_sample_inData$dbGaP_Sample_ID,]


#### getting the total number of reads 

reads <- readTotals(fdl)


write.csv(sraData, "sraData.csv")
