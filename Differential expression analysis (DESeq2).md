# Differential gene expression analysis (DESeq2)

This a general code used for gene expression analysis, the example given here is for long noncoding RNAs in luminal differentiated cells.
### 1: Admin and preparation
```
install.packages(c("RColorBrewer", "mixOmics"))
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
library("biomaRt")
ensembl    <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
                                  'chromosome_name', 'start_position', 'end_position', 
                                  'strand', 'gene_biotype'), mart = ensembl)  
 rownames(ensEMBL2id) <- ensEMBL2id$ensembl_gene_id
 raw <-MouseMammaryGland_Cleaned_rawCounts[,-1]
 rownames(raw) <- MouseMammaryGland_Cleaned_rawCounts[,1]
 annotated_raw  <- merge(raw, ensEMBL2id, by=0)
```
filter to biotypes of interest 
```
 raw_noncoding <- annotated_raw %>% dplyr::filter(annotated_raw$gene_biotype == 
                                                    "3prime_overlapping_ncRNA" |
                                                    gene_biotype == "antisense" |
                                                    gene_biotype == "bidirectional_promoter_lncRNA" |
                                                    gene_biotype == "lincRNA" |
                                                    gene_biotype == "lncRNA" |
                                                    gene_biotype == "macro_lncRNA" |
                                                    gene_biotype ==  "processed_transcript")
 
 raw_noncoding <- raw_noncoding[,c(1:453)] 
```
filter to cell type of interest
```
 sampleinfo <- MouseMammaryGland_Cleaned_MetaData
 sampleinfo_luminal <-sampleinfo %>% dplyr::filter(sampleinfo$
                                                     Cell.type =="Luminal Differentiated")
 raw_noncoding2 <- raw_noncoding %>% dplyr::select(Row.names, sampleinfo_luminal$X)
```
### 2: Running deseq2

Load data
```
 sampleinfo_luminal <- read.csv("~/Desktop/rstudio-export/sampleinfo_luminal.csv")
 rawnoncoding2 <- read.csv("~/Desktop/rstudio-export/rawnoncoding2.csv")
 library(DESeq2)
 library(ggplot2)
 metaData <- sampleinfo_luminal
 countData <- rawnoncoding2
 countData <- countData[,-c(1)]
 metaData <- metaData[,-c(1)]
 dds <- DESeqDataSetFromMatrix(countData = countData, 
                               colData = metaData,
                               design=~Stage, tidy=TRUE)
 dds <- DESeq(dds)
 res <- results(dds)
 head(results(dds,tidy=TRUE))

 summary(res)
res <- res[order(res$padj),]
head(res)
```
### 3: Plotting results
```
 par(mfrow=c(1,1))
 with(res, plot(log2FoldChange, -log10(pvalue),pch=20,main="Volcano plot",xlim=c(-3,3),ylim=c(-1,20)))
 with(subset(res,padj<0.1),points(log2FoldChange, -log10(pvalue),pch=20, col="blue"))
 with(subset(res,padj<0.1 & abs(log2FoldChange)>2),points(log2FoldChange, -log10(pvalue),pch=20,col="red"))
```
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/volcano_plot.png)

note that the names and informtation on the differentially expressed genes are saved as two separate files annotated_red.csv and annotated_blue.csv

return the names of the differentially expressed genes
```
different_blue <- subset(resdf,padj<0.1)
different_red <- subset(resdf,padj<0.1 & abs(log2FoldChange)>2)
write.csv(different_blue,file="different_blue.csv")
write.csv(different_red, file="different_red.csv")
```
annotate differential expression files
```
rownames(different_blue) <-different_blue[,1]
annotated_blue  <- merge(different_blue, ensEMBL2id, by=0)
rownames(different_red) <-different_red[,1]
annotated_red  <- merge(different_red, ensEMBL2id, by=0)
annotated_blue <- annotated_blue[,c(10,1:9,11:16)]
annotated_blue <- annotated_blue[,-c(2:3)]
annotated_red <- annotated_red[,c(10,1:9,11:16)]
annotated_red <- annotated_red[,-c(2:3)]
write.csv(annotated_blue,file="annotated_blue.csv")
write.csv(annotated_red,file="annotated_red.csv")
```
the same has been done for all genes ie not just the noncoding ones), returning two files: different_blue_all.csv and different_red_all.csv, using the same limits as before.

#### manhattan plot of differentially expressed genes using results from DESeq2

see separate file for code on generation of manhattan plot function

this shows the distribution of variation across the genome
```
 resdf <- read.csv("/data/homes/jb2220/resdf.csv")
 rownames(resdf) <- resdf[,1]
 annotated_resdf  <- merge(resdf, ensEMBL2id, by=0)
 annotated_resdf <- as.data.frame(annotated_resdf)
 annotated_resdf2 <- na.omit(annotated_resdf) #remove the NAs
 manhattan.plot(factor(annotated_resdf2$chromosome_name, levels=c(1:19, "X","Y")),main="Differential expression of noncoding RNAs in luminal differentiated > cells",
                annotated_resdf2$start_position,annotated_resdf2$pvalue,sig.level=5e-3)
```
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-11%20at%2010.17.25.png)

