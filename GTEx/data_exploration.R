library(ngsReports)
library(dplyr)
library(magrittr)
library(stringr)
library(kableExtra)
library(AMR)
library(summarytools)
library(pander)


### Load all the data 
gtex_data<- read.table("Counts_matrix/featurecounts_GTEX_Rmatrix.txt", header=TRUE, row.names = 1)
inHouse_edata <- read.table("Counts_matrix/featurecounts_Rmatrix.txt", header=TRUE, row.names = 1)
pdata_gtex <- read.csv("Pdata_GTEX.csv", header=TRUE)
pdata_all <- read.csv("pdata_all.csv", header=TRUE)


### using only RNA-seq files 

gtex_data <- t(gtex_data)
gtex_edata <- gtex_data[rownames(gtex_data) %in% pdata_gtex$Run,]
gtex_edata <- t(gtex_edata)
dim(gtex_edata)


pdata_gtex <- pdata_gtex[order(match(rownames(pdata_gtex),colnames(gtex_edata))),]



#### change column names 
colnames(inHouse_edata) <- c("UPF3B_03", "UPF2_01", "UPF2_02", "UPF2_03", "UPF2_04", "Control_01", "Control_02", 
                             "Control_03", "Control_04", "Control_05", "Control_06", "Control_07", "UPF3B_01", "UPF3B_02")


#### merge all files together 
combined <- merge(as.data.frame(gtex_edata), as.data.frame(inHouse_edata), by = "row.names", sort=FALSE)

combined <- combined %>%
  column_to_rownames("Row.names")
combined

dim(combined)

pdata_all <- pdata_all %>%
  column_to_rownames('Run')

genes2keep <- cpm(combined) %>%
  is_greater_than(1) %>%
  rowSums() %>%
  is_greater_than(3)
combined <- combined[genes2keep,]

pdata_all <- pdata_all[order(match(rownames(pdata_all),colnames(combined))),]

pca <-prcomp(t(cpm(combined, log = T)))
pcaPlot <- pca$x[,1:5] %>%
  as.data.frame() %>%
  rownames_to_column('Run') %>%
  cbind(pdata_all[.$Run,]) %>%
  ggplot(aes(PC4, PC5, key = Run, shape = Condition)) + geom_point(size = 11.5) + theme_bw() + 
geom_text(aes(label=rownames(pdata_all), color = "black"))
pcaPlot

dim(combined)


#### get rid of the weird gtex sample


combined_dropped <- subset(combined, select= -c(SRR820914))

pdata_all <- pdata_all[rownames(pdata_all) %in% colnames(combined_dropped),]
pdata_all <- pdata_all[order(match(rownames(pdata_all),colnames(combined_dropped))),]

pca <-prcomp(t(cpm(combined_dropped, log = T)))
pcaPlot <- pca$x[,1:5] %>%
  as.data.frame() %>%
  rownames_to_column('Run') %>%
  cbind(pdata_all[.$Run,]) %>%
  ggplot(aes(PC3, PC4, key = Run, shape = Condition)) + geom_point(size = 11.5) + theme_bw() + 
  geom_text(aes(label=rownames(pdata_all), color = "black"))
pcaPlot



dds <- DESeqDataSetFromMatrix(countData = combined_dropped,
                              colData = pdata_all,
                              design= ~Condition + SEX + AGE + BATCH)
dds = estimateSizeFactors(object = dds)

dds <- DESeq(dds)

resultsNames(dds)

Upf2_vs_Controls <- DESeq2::results(dds, contrast = c("Condition", "UPF2.", "Control"), pAdjustMethod = "BH")
Upf2_vs_Controls <- Upf2_vs_Controls[order(Upf2_vs_Controls$padj), ]
summary(Upf2_vs_Controls, alpha = 0.05)
Upf2_vs_Controls$symbol = mapIds(org.Hs.eg.db,
                                 keys=row.names(Upf2_vs_Controls), 
                                 column="SYMBOL",
                                 keytype="ENSEMBL",
                                 multiVals="first")


Upf2_vs_Controls <- Upf2_vs_Controls %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  dplyr::filter(padj < 0.05) %>%
  column_to_rownames('gene')

head(Upf2_vs_Controls, 25)

Upf2P_vs_Controls <- DESeq2::results(dds, contrast = c("Condition", "UPF2P", "Control"), pAdjustMethod = "BH")
Upf2P_vs_Controls <- Upf2P_vs_Controls[order(Upf2P_vs_Controls$padj), ]
summary(Upf2P_vs_Controls, alpha = 0.05)


head(Upf2P_vs_Controls, 30)

Upf2P_vs_Controls$symbol = mapIds(org.Hs.eg.db,
                                  keys=row.names(Upf2P_vs_Controls), 
                                  column="SYMBOL",
                                  keytype="ENSEMBL",
                                  multiVals="first")

head(Upf2_vs_Controls, 30)


Upf3B_vs_Controls <- DESeq2::results(dds, contrast = c("Condition", "UPF3B", "Control"), pAdjustMethod = "BH")
Upf3B_vs_Controls <- Upf3B_vs_Controls[order(Upf3B_vs_Controls$padj), ]
summary(Upf3B_vs_Controls, alpha = 0.05)


head(Upf3B_vs_Controls, 30)

Upf3B_vs_Controls$symbol = mapIds(org.Hs.eg.db,
                                  keys=row.names(Upf3B_vs_Controls), 
                                  column="SYMBOL",
                                  keytype="ENSEMBL",
                                  multiVals="first")


Upf3B_vs_Controls <- Upf3B_vs_Controls %>%
  as.data.frame() %>%
  rownames_to_column('gene') %>%
  dplyr::filter(padj < 0.05) %>%
  column_to_rownames('gene')
