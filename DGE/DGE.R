library(tibble)
library(magrittr)
library(limma)
library(edgeR)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(org.Hs.eg.db)
library(VennDiagram)
#------------------------ Loading all the data-------------------#

counts_data <- read.table("files/featurecounts_Rmatrix.txt", header= TRUE, row.names = 1)
head(counts_data)

colnames(counts_data) <- c("16-PX1hN", 
                           "17-PX2hN", "18-PX5hN", "19-3AG1C3hN","20-3AG2C5hN",
                           "21-3AG3C14hN","22-3BG1C2hN", "23-3BG2C8hN","24_HI","25_HI",  "26-PX1hE", 
                           "27-PX2hE", "28-PX3hE", "29-PX5hE",
                           "30-3AG1C3hE", "31-3AG2C5hE", 
                           "32-3AG3C14hE", "33-3BG1C2hE", "34-3BG2C8hE",
                           "35-3BG3C12hE")

counts_ES <- counts_data %>%
  dplyr::select("26-PX1hE", 
                "27-PX2hE", "28-PX3hE", "29-PX5hE",
                "30-3AG1C3hE", "31-3AG2C5hE", 
                "32-3AG3C14hE", "33-3BG1C2hE", "34-3BG2C8hE",
                "35-3BG3C12hE")
head(counts_ES)


pdata <- read.csv("~/Urwah/Debrah/Salmon/ES_cells/labels.txt", header= TRUE)

pdata <- pdata %>% 
  column_to_rownames("Sample")

rownames(counts_ES) <- gsub("\\_[0-9]*$", "", rownames(counts_ES))
rownames(counts_ES) <- gsub("\\.[0-9]*$", "", rownames(counts_ES))

head(counts_ES)
#-----------------------------------------------------------------#
#-----------------------Exploratory analysis----------------------#
#-----------------------------------------------------------------#

y <- DGEList(counts_ES)

### Plotting 
plotMDS(y)

pca <-prcomp(t(cpm(y$counts, log = T)))
pca$x[,1:4] %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  cbind(pdata[.$Sample,]) %>%
  ggplot(aes(PC1, PC2, key = Sample, color = Condition, shape = CollectionDate)) + geom_point(size = 11.5) + theme_bw()


genes2keep <- cpm(y$counts) %>%
  is_greater_than(1) %>%
  rowSums() %>%
  is_greater_than(3)
y <- y[genes2keep,]

dim(y$counts)
pca <-prcomp(t(cpm(y$counts, log = T)))
pca$x[,1:4] %>%
  as.data.frame() %>%
  rownames_to_column('Sample') %>%
  cbind(pdata[.$Sample,]) %>%
  ggplot(aes(PC1, PC2, key = Sample, color = Condition, shape = CollectionDate)) + geom_point(size = 11.5) + theme_bw()

plotMDS(y)


y <- calcNormFactors(y)


### DE analysis 

design <- model.matrix(~Condition + CollectionDate, data =pdata)
design
y <- voom(y, design, plot = T)

fit <- lmFit(y, design)
head(coef(fit))
colnames(coef(fit))

rm(tmp)
tmp <- eBayes(fit)
top.table <- topTable(tmp, sort.by = "P", n = Inf, coef = "ConditionUpf3bKO")

top.table
summary(decideTests(top.table))


?topTable
rownames(top.table) <- gsub("\\.[0-9]*$","", rownames(top.table))
head(top.table)


top.table$symbol = mapIds(org.Hs.eg.db,
                                  keys=row.names(top.table), 
                                  column="SYMBOL",
                                  keytype="ENSEMBL",
                                  multiVals="first")
hist(top.table$P.Value)
head(top.table,100)

sig_genes <- top.table %>%
  dplyr::filter(adj.P.Val < 0.05)

dim(sig_genes)

rm(y)
#---- Voom with quality weights 
v <- voom(y, design = design)
plotMDS(v)
vfit <- lmFit(v)
vfit <- eBayes(vfit)
topTable(vfit, coef = "ConditionUpf3bKO", sort.by = "P")
top <- topTable(vfit,coef=3,number=Inf,sort.by="P")
top
sum(top$adj.P.Val<0.05)

rownames(top) <- gsub("\\.[0-9]*$","", rownames(top))
head(top)


top$entrez = mapIds(org.Hs.eg.db,
                                  keys=row.names(top), 
                                  column="ENTREZID",
                                  keytype="ENSEMBL",
                                  multiVals="first")
top$symbol = mapIds(org.Hs.eg.db,
                                  keys=row.names(top), 
                                  column="SYMBOL",
                                  keytype="ENSEMBL",
                                  multiVals="first")

top
################################################################################
################################################################################
vwts <- voomWithQualityWeights(y, design=design, normalization="none", plot=TRUE)
vfit2 <- lmFit(vwts)
vfit2 <- eBayes(vfit2)
top2<- topTable(vfit2,coef=3,sort.by="P")
sum(top2$adj.P.Val<0.05)
top2

rownames(top2) <- gsub("\\.[0-9]*$","", rownames(top2))
head(top2)


top2$symbol = mapIds(org.Hs.eg.db,
                          keys=row.names(top2), 
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
top2

##### Gene set testing 

go <- goana(vfit, species = "Hs")
?goana


####################################################################################
####################################################################################
## EdgeR analysis 
##

y <- DGEList(counts_ES)


genes2keep <- cpm(y$counts) %>%
  is_greater_than(1) %>%
  rowSums() %>%
  is_greater_than(3)
y <- y[genes2keep,]

barplot(y$samples$lib.size,names=colnames(y),las=2)



dgList <- estimateGLMCommonDisp(y, design=design)
dgList <- estimateGLMTagwiseDisp(dgList, design=design)
plotBCV(dgList)

fit <- glmFit(dgList, design)

lrt <- glmLRT(fit, coef = 3)
summary(decideTests(lrt))

######################### edgeR quasi F test 
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)
y$common.dispersion

fit <- glmQLFit(y, design, robust=TRUE)
fit

### UPF3B KO
Upf3bKO <- glmQLFTest(fit, coef = 3)
summary(decideTests(Upf3bKO))

####
Upf3bKO_treat <- glmTreat(fit, coef = 3)
summary(decideTests(Upf3bKO_treat))
Upf3bKO_treat_table <- topTags(Upf3bKO_treat)

###############
### annotate 

Upf3bKO_treat_table$symbol = mapIds(org.Hs.eg.db,
                          keys=row.names(Upf3bKO_treat_table), 
                          column="SYMBOL",
                          keytype="ENSEMBL",
                          multiVals="first")
Upf3bKO_treat_table$test
