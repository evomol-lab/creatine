library(GEOquery)
library(dplyr)
library(ggplot2)
library(ggrepel)
## change my_id to be the dataset that you want.
my_id <- "GSE6344"
gse <- getGEO(my_id)
length(gse)
gse1 <- gse[[1]]
gse2 <- gse[[2]]
exprs(gse1) <- log2(exprs(gse1))
exprs(gse2) <- log2(exprs(gse2))
boxplot(exprs(gse1),outline=FALSE)
boxplot(exprs(gse2),outline=FALSE)
#library(dplyr)
sampleInfo1 <- pData(gse1)
sampleInfo2 <- pData(gse2)
sampleInfo1 <- select(sampleInfo1, "title","characteristics_ch1")
sampleInfo2 <- select(sampleInfo2, "title","characteristics_ch1")
sampleInfo2 <- rename(sampleInfo2, "Sample"="title", "State"="characteristics_ch1")
sampleInfo1 <- rename(sampleInfo1, "Sample"="title", "State"="characteristics_ch1")
View(sampleInfo1)
View(sampleInfo2)
sampleInfo1$Sample <- sub(",.*", "", sampleInfo1$Sample)
sampleInfo1$State <- sub(",.*", "", sampleInfo1$State)
sampleInfo2$Sample <- sub(",.*", "", sampleInfo2$Sample)
sampleInfo2$State <- sub(",.*", "", sampleInfo2$State)
library(pheatmap)
corMatrix1 <- cor(exprs(gse1),use="c")
pheatmap(corMatrix1, annotation_col=sampleInfo1)
corMatrix2 <- cor(exprs(gse2),use="c")
pheatmap(corMatrix2, annotation_col=sampleInfo2)
#library(ggplot2)
#library(ggrepel)
pca1 <- prcomp(t(exprs(gse1)))
pca2 <- prcomp(t(exprs(gse2)))
## Join the PCs to the sample information
cbind(sampleInfo1, pca1$x) %>%
    ggplot(aes(x = PC1, y=PC2, col=Sample,label=paste("", State))) + geom_point() + geom_text_repel()
cbind(sampleInfo2, pca2$x) %>%
    ggplot(aes(x = PC1, y=PC2, col=Sample,label=paste("", State))) + geom_point() + geom_text_repel()
# DEs for 4 designs
library(limma)
# Based on the column Sample
#design1 <- model.matrix(~0+sampleInfo1$Sample)
#design2 <- model.matrix(~0+sampleInfo2$Sample)
#design1
#design2
# Based on the column State
design21 <- model.matrix(~0+sampleInfo1$State)
design22 <- model.matrix(~0+sampleInfo2$State)
design21
design22
## the column names are a bit ugly, so we will rename
#colnames(design1) <- c("Stage1","Stage2", "Tumor")
#colnames(design2) <- c("Stage1","Stage2")
colnames(design21) <- c("Normal","Tumor")
colnames(design22) <- c("Normal","Tumor")

summary(exprs(gse1))
summary(exprs(gse2))

## calculate median expression level
cutoff1 <- median(exprs(gse1))
cutoff2 <- median(exprs(gse2))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed1 <- exprs(gse1) > cutoff1
is_expressed2 <- exprs(gse2) > cutoff2
## Identify genes expressed in more than 2 samples

keep1 <- rowSums(is_expressed1) > 2
keep2 <- rowSums(is_expressed2) > 2

## check how many genes are removed / retained.
table(keep1)
table(keep2)

## subset to just those expressed genes
g1se <- gse1[keep1,]
g2se <- gse2[keep2,]

## Fit to all Stages/States
#fit1 <- lmFit(exprs(g1se), design1)
#head(fit1$coefficients)
#fit2 <- lmFit(exprs(g2se), design2)
#head(fit2$coefficients)
fit21 <- lmFit(exprs(g1se), design21)
head(fit21$coefficients)
fit22 <- lmFit(exprs(g2se), design22)
head(fit22$coefficients)

## Now the contrasts
#contrasts1 <- makeContrasts(Stage2 - Stage1, levels=design1)
#contrasts2 <- makeContrasts(Stage2 - Stage1, levels=design2)
contrasts21 <- makeContrasts(Tumor - Normal, levels=design21)
contrasts22 <- makeContrasts(Tumor - Normal, levels=design22)

#fit1 <- contrasts.fit(fit1, contrasts1)
#fit2 <- contrasts.fit(fit2, contrasts2)
fit21 <- contrasts.fit(fit21, contrasts21)
fit22 <- contrasts.fit(fit22, contrasts22)

#fit1 <- eBayes(fit1)
#fit2 <- eBayes(fit2)
fit21 <- eBayes(fit21)
fit22 <- eBayes(fit22)

#topTable(fit1)
#topTable(fit2)
topTable(fit21)
topTable(fit22)

#table(decideTests(fit1))
#table(decideTests(fit2))
table(decideTests(fit21))
table(decideTests(fit22))

# Visualization of DE Results

## Getting Annotations
anno1 <- fData(gse1)
anno2 <- fData(gse2)

anno1 <- select(anno1,"GB_ACC", "Gene Symbol")
anno2 <- select(anno2,"GB_ACC", "Gene Symbol")

#fit1$genes <- anno1
#topTable(fit21)
#fit2$genes <- anno2
#topTable(fit21)

fit21$genes <- anno1
topTable(fit21)
fit22$genes <- anno2
topTable(fit21)


full_results21 <- topTable(fit21, number=Inf)
full_results21 <- tibble::rownames_to_column(full_results21,"ID")
full_results22 <- topTable(fit22, number=Inf)
full_results22 <- tibble::rownames_to_column(full_results22,"ID")

library(ggplot2)
ggplot(full_results21,aes(x = logFC, y=B)) + geom_point()
ggplot(full_results22,aes(x = logFC, y=B)) + geom_point()
p_cutoff <- 0.05
fc_cutoff <- 1
full_results21 %>% 
    mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
    ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()
full_results22 %>% 
    mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
    ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()
library(ggrepel)
p_cutoff <- 0.05
fc_cutoff <- 1
topN <- 20
full_results21 %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, Gene.Symbol,"")) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")
full_results22 %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, Gene.Symbol,"")) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")

fullresultsfiltered21 <- filter(full_results21, adj.P.Val < 0.05, abs(logFC) > 1)
fullresultsfiltered22 <- filter(full_results22, adj.P.Val < 0.05, abs(logFC) > 1)
genes_of_interest <- c('SLC2A4', 'IGF1', 'GATM', 'SLC6A8', 'AKT1', 'AKT2', 'AKT3', 'PSMB5', 'PSMD3', 'SGK1', 'PRPS1', 'CKB', 'CKM', 'CKMT1A', 'CKMT1B', 'CKMT2', 'SLC16A12')

ids_of_interest21 <-  filter(full_results21,Gene.Symbol %in% genes_of_interest) %>% 
    pull(ID)
gene_names21 <-  filter(full_results21,Gene.Symbol %in% genes_of_interest) %>% 
    pull(Gene.Symbol)
ids_of_interest22 <-  filter(full_results22,Gene.Symbol %in% genes_of_interest) %>% 
    pull(ID)
gene_names22 <-  filter(full_results22,Gene.Symbol %in% genes_of_interest) %>% 
    pull(Gene.Symbol)
gene_matrix21 <- exprs(gse1)[ids_of_interest21,]
gene_matrix22 <- exprs(gse2)[ids_of_interest22,]
pheatmap(gene_matrix21,
         labels_row = gene_names21,
         scale="row", annotation_col=sampleInfo1)
pheatmap(gene_matrix22,
         labels_row = gene_names22,
         scale="row", annotation_col=sampleInfo2)
ids_of_interestfiltered21 <-  filter(fullresultsfiltered21,Gene.Symbol %in% genes_of_interest) %>% 
    pull(ID)
gene_namesfiltered21 <-  filter(fullresultsfiltered21,Gene.Symbol %in% genes_of_interest) %>% 
    pull(Gene.Symbol)
ids_of_interestfiltered22 <-  filter(fullresultsfiltered22,Gene.Symbol %in% genes_of_interest) %>% 
    pull(ID)
gene_namesfiltered22 <-  filter(fullresultsfiltered22,Gene.Symbol %in% genes_of_interest) %>% 
    pull(Gene.Symbol)
gene_matrixfiltered21 <- exprs(gse1)[ids_of_interestfiltered21,]
gene_matrixfiltered22 <- exprs(gse2)[ids_of_interestfiltered22,]
pheatmap(gene_matrixfiltered21,
         labels_row = gene_namesfiltered21,
         scale="row", annotation_col=sampleInfo1)
pheatmap(gene_matrixfiltered22,
         labels_row = gene_namesfiltered22,
         scale="row", annotation_col=sampleInfo2)