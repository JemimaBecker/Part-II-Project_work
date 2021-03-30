### Prep

load libraries 
```
> install.packages(c("RColorBrewer", "mixOmics"))
> source("http://bioconductor.org/biocLite.R")
> biocLite("edgeR")

> if(!requireNamespace("BiocManager"))
>   install.packages("BiocManager")
> BiocManager::install(c("limma", "edgeR","Glimma", "org.Mm.eg.db", "gplots", "RColorBrewer", "NMF", "BiasedUrn"))

> library(edgeR)
> library(limma)
> library(Glimma)
> library(org.Mm.eg.db)
> library(gplots)
> library(RColorBrewer)
> library(NMF)
```
### create object that is just the counts data for lncRNAS

1: loading ensembl annotations
```
> MouseMammaryGland_Cleaned_rawCounts <- read.csv("/data/homes/rsh46/MammaryGlandData/MouseMammaryGland_Cleaned_rawCounts.csv")
> library("biomaRt")

> message("+-------------------------------------------------------------------------------")
> message("+ Use ensEMBL Annotations")
> message("+-------------------------------------------------------------------------------")

> ensembl    <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
> ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
>                                  'chromosome_name', 'start_position', 'end_position', 
>                                  'strand', 'gene_biotype'), mart = ensembl)  
> rownames(ensEMBL2id) <- ensEMBL2id$ensembl_gene_id

> message("+-------------------------------------------------------------------------------")
> message("+ Merge Expression and ensEMBL Annotations")
> message("+-------------------------------------------------------------------------------")
```
2: adding annotations to raw counts and filtering 
```
> raw <-MouseMammaryGland_Cleaned_rawCounts[,-1]
> rownames(raw) <- MouseMammaryGland_Cleaned_rawCounts[,1]
> annotated_raw  <- merge(FPKM_unG, ensEMBL2id, by=0)

> raw_noncoding <- annotated_raw %>% dplyr::filter(annotated_raw$gene_biotype == 
>                                                    "3prime_overlapping_ncRNA" |
>                                                    gene_biotype == "antisense" |
>                                                    gene_biotype == "bidirectional_promoter_lncRNA" |
>                                                    gene_biotype == "lincRNA" |
>                                                    gene_biotype == "lncRNA" |
>                                                    gene_biotype == "macro_lncRNA" |
>                                                    gene_biotype ==  "processed_transcript")


> raw_noncoding <- raw_noncoding[,c(1:453)] #remove other annotations
```
### prepare and organise data
```
> sampleinfo <- MouseMammaryGland_Cleaned_MetaData
> seqdata <- raw_noncoding
> countdata <- seqdata[,-c(1)]
> rownames(countdata) <- seqdata[,1]
> table(colnames(countdata)==sampleinfo$X)
```
convert counts to DGEList object 
```
> y <- DGEList(countdata)
> group <- paste(sampleinfo$Cell.type, sampleinfo$Stage, sep=".")
> group <- factor(group)
> y$samples$group <- group
```
filtering lowly expressed genes
```
> myCPM <- cpm(countdata)
> thresh <- myCPM > 0.5
> keep <- rowSums(thresh) > 2
> plot(myCPM[,1],countdata[,1])
```
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-04%20at%2021.39.11.png)

### multidimensional scaling plots
```
> plotMDS(y)
> par(mfrow=c(1,2))
> col.cell <- c("purple","orange", "red", "blue", "yellow", "green","pink")[sampleinfo$Cell.type]
> data.frame(sampleinfo$Cell.type,col.cell)

> plotMDS(y,dim=c(3,4),col=col.cell,pch=16)
> legend("topleft",fill=c("purple","orange", "red", "blue", "yellow", "green","pink"),legend=levels(sampleinfo$Cell.type))
> title("Cell type")
```
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-04%20at%2021.38.57.png)
