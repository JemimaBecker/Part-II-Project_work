### 1: Admin and preparation

> install.packages(c("RColorBrewer", "mixOmics"))
> source("http://bioconductor.org/biocLite.R")
> biocLite("edgeR")
> library("biomaRt")
> ensembl    <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
> ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
>                                  'chromosome_name', 'start_position', 'end_position', 
>                                  'strand', 'gene_biotype'), mart = ensembl)  
> rownames(ensEMBL2id) <- ensEMBL2id$ensembl_gene_id
> raw <-MouseMammaryGland_Cleaned_rawCounts[,-1]
> rownames(raw) <- MouseMammaryGland_Cleaned_rawCounts[,1]
> annotated_raw  <- merge(raw, ensEMBL2id, by=0)

filter to biotypes of interest 

> raw_noncoding <- annotated_raw %>% dplyr::filter(annotated_raw$gene_biotype == 
>                                                    "3prime_overlapping_ncRNA" |
>                                                    gene_biotype == "antisense" |
>                                                    gene_biotype == "bidirectional_promoter_lncRNA" |
>                                                    gene_biotype == "lincRNA" |
>                                                    gene_biotype == "lncRNA" |
>                                                    gene_biotype == "macro_lncRNA" |
>                                                    gene_biotype ==  "processed_transcript")
> 
> raw_noncoding <- raw_noncoding[,c(1:453)] 

filter to cell type of interest

> sampleinfo <- MouseMammaryGland_Cleaned_MetaData
> sampleinfo_luminal <-sampleinfo %>% dplyr::filter(sampleinfo$
>                                                     Cell.type =="Luminal Differentiated")
> raw_noncoding2 <- raw_noncoding %>% dplyr::select(Row.names, sampleinfo_luminal$X)

### 2: doing deseq2

load data

> sampleinfo_luminal <- read.csv("~/Desktop/rstudio-export/sampleinfo_luminal.csv")
> rawnoncoding2 <- read.csv("~/Desktop/rstudio-export/rawnoncoding2.csv")
> library(DESeq2)
> library(ggplot2)
> metaData <- sampleinfo_luminal
> countData <- rawnoncoding2
> countData <- countData[,-c(1)]
> metaData <- metaData[,-c(1)]
> dds <- DESeqDataSetFromMatrix(countData = countData, 
>                               colData = metaData,
>                               design=~Stage, tidy=TRUE)
> dds <- DESeq(dds)
> res <- results(dds)
> head(results(dds,tidy=TRUE))

> summary(res)

out of 5231 with nonzero total read count

adjusted p-value < 0.1

LFC > 0 (up)       : 163, 3.1%

LFC < 0 (down)     : 244, 4.7%

outliers [1]       : 0, 0%

low counts [2]     : 2937, 56%

(mean count < 1)

[1] see 'cooksCutoff' argument of ?results

[2] see 'independentFiltering' argument of ?results

> res <- res[order(res$padj),]
> head(res)

### 3: Plotting results

> par(mfrow=c(1,1))
> with(res, plot(log2FoldChange, -log10(pvalue),pch=20,main="Volcano plot",xlim=c(-3,3),ylim=c(-1,20)))
> with(subset(res,padj<0.1),points(log2FoldChange, -log10(pvalue),pch=20, col="blue"))
> with(subset(res,padj<0.1 & abs(log2FoldChange)>2),points(log2FoldChange, -log10(pvalue),pch=20,col="red"))

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/volcano_plot.png)


return the names of the differentially expressed genes
> different_blue <- subset(resdf,padj<0.1)
> different_red <- subset(resdf,padj<0.1 & abs(log2FoldChange)>2)
> write.csv(different_blue,file="different_blue.csv")
> write.csv(different_red, file="different_red.csv")

annotate differential expression files

> rownames(different_blue) <-different_blue[,1]
> annotated_blue  <- merge(different_blue, ensEMBL2id, by=0)
> rownames(different_red) <-different_red[,1]
> annotated_red  <- merge(different_red, ensEMBL2id, by=0)
> annotated_blue <- annotated_blue[,c(10,1:9,11:16)]
> annotated_blue <- annotated_blue[,-c(2:3)]
> annotated_red <- annotated_red[,c(10,1:9,11:16)]
> annotated_red <- annotated_red[,-c(2:3)]
> write.csv(annotated_blue,file="annotated_blue.csv")
> write.csv(annotated_red,file="annotated_red.csv")
