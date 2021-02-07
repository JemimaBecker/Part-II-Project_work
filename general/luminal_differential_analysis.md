## 1: Preparation
### 1: installing packages

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

### 2: creating object that is raw counts data for noncoding genes in luminal differentiated and progenitor cells

> load ensembl annotations

> library("biomaRt")
> ensembl    <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
> ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
>                                  'chromosome_name', 'start_position', 'end_position', 
>                                  'strand', 'gene_biotype'), mart = ensembl)  
> rownames(ensEMBL2id) <- ensEMBL2id$ensembl_gene_id
> raw <-MouseMammaryGland_Cleaned_rawCounts[,-1]
> rownames(raw) <- MouseMammaryGland_Cleaned_rawCounts[,1]
> annotated_raw  <- merge(FPKM_unG, ensEMBL2id, by=0)

filter to specific gene biotypes

> raw_noncoding <- annotated_raw %>% dplyr::filter(annotated_raw$gene_biotype == 
>                                                    "3prime_overlapping_ncRNA" |
>                                                    gene_biotype == "antisense" |
>                                                    gene_biotype == "bidirectional_promoter_lncRNA" |
>                                                    gene_biotype == "lincRNA" |
>                                                    gene_biotype == "lncRNA" |
>                                                    gene_biotype == "macro_lncRNA" |
>                                                    gene_biotype ==  "processed_transcript")


> raw_noncoding <- raw_noncoding[,c(1:453)]

filter to cell types of interest

> sampleinfo <- MouseMammaryGland_Cleaned_MetaData
> sampleinfo_luminal <-sampleinfo %>% dplyr::filter(sampleinfo$
>                                                     Cell.type =="Luminal Differentiated" | 
>                                                     Cell.type == "Luminal Progenitors")

> raw_noncoding2 <- raw_noncoding %>% dplyr::select(Row.names, sampleinfo_luminal$X)
> seqdata <- raw_noncoding2
> countdata <- seqdata [,-c(1)]
> rownames(countdata) <- seqdata[,1]
> table> (colnames(countdata)==sampleinfo_luminal$X> )

#### 3: convert counts to DGEList object

> y <- DGEList(countdat> a)
> group <- paste(sample> info_luminal$Cell.type, sampleinfo_luminal$Stage, sep=".")
> group <- factor(grou> p)
> y$samples$group <- g> rou> > > p

## 2: Hierarchical clustering

> var_genes <- apply(logcounts,1,var)
> head(var_genes)

return list of top 100 genes that show the greatest variance

> select_var <- names(sort(var_genes, decreasing = TRUE))[1:100]
> selectvar_luminal <- as.data.frame(select_var)

annotate this list

> ensembl    <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
> ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
>                                  'chromosome_name', 'start_position', 'end_position', 
>                                  'strand', 'gene_biotype'), mart = ensembl)  
> rownames(ensEMBL2id) <- ensEMBL2id$ensembl_gene_id
> rownames(selectvar_luminal) <- selectvar_luminal[,1]
> selectvar_luminal_annotated  <- merge(selectvar_luminal, ensEMBL2id, by=0)

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-07%20at%2015.58.13.png)

a problem here: some of the gene biotypes can apply to both noncoding and protien coding RNAs, so by filtering by biotype I have to choose between including ALL ncRNAs (and including a few protein coding genes) or only including ncRNAs and missing a few of them. >:( shit 

> highly_variable_lcpm <- logcounts[select_var,]
> dim(highly_variable_lcpm)

plotting 

> mypalette <- brewer.pal(11,"RdYlBu")
> morecols <- colorRampPalette(mypalette)
> col.cell <- c("purple","pink","orange","red","blue","green","yellow")[sampleinfo_luminal$Cell.type]

had an error trying to run col.cell - will try again at some point

> heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none",main="Top 100 most variable genes across Luminal samples",
>           scale="row",cex.lab=0.6)
> legend(inset=c(-0.2,0),cex=0.5,"bottomleft",fill=c("purple","pink","orange","red","blue","green","yellow"),legend=levels(sampleinfo$Cell.type))

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-07%20at%2015.56.05.png)
