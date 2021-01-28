> install.packages("htmltools")
> library(htmltools)
> #source("https://bioconductor.org/biocLite.R")
> biocLite("DESeq2")

> library( "DESeq2" )
> library(ggplot2)
> library("BiocParallel")
> register(MulticoreParam(4))
> countData <- MouseMammaryGland_Cleaned_rawCounts
> head(countData)
> metaData <- MouseMammaryGland_Cleaned_MetaData
> dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~Stage, tidy = TRUE)
> dds <- DESeq(dds)
> ?DESeq
> res <- results(dds)
> head(results(dds, tidy=TRUE))
> summary(res)
> res <- res[order(res$padj),]
> head(res)
