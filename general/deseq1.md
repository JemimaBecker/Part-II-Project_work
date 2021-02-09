# Running DeSeq2
29-01-2021 Downloaded data and ran on desktop R studio after overwhelming jeremy
> #install.packages("htmltools")
> #library(htmltools)
> #source("https://bioconductor.org/biocLite.R")
> #biocLite("DESeq2")

> #admin stuff
> library( "DESeq2" )
> library(ggplot2)
> library("BiocParallel")
> register(MulticoreParam(4))

# by stage
> countData <- Home_rawcounts
> head(countData)
> metaData <- metadata.Name
> dds <- DESeqDataSetFromMatrix(countData=countData, 
>                               colData=metaData, 
>                               design=~Stage, tidy = TRUE)


> dds <- DESeq(dds, parallel=T)
> res <- results(dds, parallel=T)
> head(results(dds, tidy=TRUE))
> summary(res)
> res <- res[order(res$padj),]
> head(res)

> resSig <- res[which(res$padj < 0.1),]
> head(resSig[order(resSig$log2FoldChange),])

# Plotting and visualisation
plotCounts(dds, gene ="ENSMUSG00000000031", intgroup = "Stage",main="H19 by stage")
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-01-30%20at%2015.08.06.png)
> with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="volcano plot", xlim=c(-3,3)))
> with(subset(res,padj<0.1), points(log2FoldChange, -log10(pvalue),pch=20, col="blue"))
> with(subset(res, padj<0.1 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-01-30%20at%2015.05.15.png)
> vsdata <- vst(dds, blind=FALSE)
> plotPCA(vsdata, intgroup="Stage")
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-01-28%20at%2014.58.56.png)
