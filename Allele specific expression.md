# Allele specific expression

Code kindly adapted from Martin Limback-Stokin
```
load.Rdata(file="/home/rsh46/MammaryGlandData/ASE/MouseMammaryGland_ASE_NormalisedReadCounts_RAW.Rdata","ASE")
ensembl    <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
                                 'chromosome_name', 'start_position', 'end_position', 
                                 'strand', 'gene_biotype'), mart = ensembl)  
rownames(ensEMBL2id) <- ensEMBL2id$ensembl_gene_id
ASE_annotated  <- merge(ASE, ensEMBL2id, by=0)
lncRNAs_ASE <- ASE_annotated %>% dplyr::filter(ASE_annotated$gene_biotype ==
                                                                   "3prime_overlapping_ncRNA" |
                                                                   gene_biotype == "antisense" |
                                                                   gene_biotype == "bidirectional_promoter_lncRNA" |
                                                                   gene_biotype == "lincRNA" |
                                                                   gene_biotype == "lncRNA" |
                                                                   gene_biotype == "macro_lncRNA" |
                                                                   gene_biotype ==  "processed_transcript")
#need to split by cell type
luminal_lncRNAs_ASE <- lncRNAs_ASE %>% dplyr::select(Row.names,wk8_meta_data_luminal$X)
wk8_metaData_luminal <-MouseMammaryGland_Cleaned_MetaData %>% dplyr::filter(MouseMammaryGland_Cleaned_MetaData$
                                                    Cell.type =="Luminal Differentiated")
```
## Allelic expression & Visualisation
```
library(edgeR)
library(limma)
library(Glimma)
library(tidyverse)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)
library(biomaRt)
library(reshape)
library("ggpubr")
ase_counts <- load("/home/rsh46/MammaryGlandData/ASE/MouseMammaryGland_ASE_NormalisedReadCounts_RAW.Rdata")
colnames(normCounts)
baseDir <- "/home/rsh46/MammaryGlandData/"
metadata  <- read.delim(paste0(baseDir, "MouseMammaryGland_Cleaned_MetaData.tsv"), header=TRUE)
LD_meta  <- subset(metadata, Cell.type=="Luminal Differentiated")
```

Format metadata
```
LD_meta1 <- as.data.frame(LD_meta)
LD_meta2 <- as.data.frame(LD_meta)
LD_meta1$sampleNames <- paste0(LD_meta1$sampleNames, ".", "g1")
LD_meta2$sampleNames <- paste0(LD_meta2$sampleNames, ".", "g2")
full_LD_meta <- rbind(LD_meta1, LD_meta2)
rownames(full_LD_meta) <- full_LD_meta$sampleNames
samples <- rownames(LD_meta)
samples1 <- paste0(samples, ".", "g1")
samples2 <- paste0(samples, ".", "g2")
all_samples <- c(samples1, samples2) 
LD_normCounts <- normCounts[, all_samples]
dim(LD_normCounts)
countdata <- LD_normCounts
all(colnames(countdata)==full_LD_meta$sampleNames)
```
Filter to remove lowly expressed genes
Use CPM function from edgeR to generate CPM values and filter. Converting to CPMs normalises for different sequencing depths of each sample.
```
myCPM <- cpm(countdata)
head(myCPM)
```
Find values in "myCPM" greater than 0.5 (can be changed) This produces a logical matrix with TRUEs and FALSEs.
```
thresh <- myCPM > 0.5
table(rowSums(thresh))
```
Want to keep genes that have at least 4 TRUES in each row of thresh (we have 8 biological replicates)
```
keep <- rowSums(thresh) >= 4
summary(keep)
counts.keep <- countdata[keep,]
dim(countdata)
dim(counts.keep)
```
Convert counts to `DGEList` object, which will also convert to fpkm and include metada 
read in gene lengths (last step subsets to only the genes that we have kept in "counts.keep")
```
edgeR_gene_lengths <- read.table(paste0(baseDir, "SLX-18042.D701rna_D505rna.H3YKJDSXY.s_4.r_1.sorted_gene.featureCounts.txt"), header=TRUE, row.names=1)
edgeR_gene_lengths <- edgeR_gene_lengths[ order(row.names(edgeR_gene_lengths)), c("Length","gene_name")]
edgeR_gene_lengths <- edgeR_gene_lengths[rownames(counts.keep),]
full_LD_meta <- full_LD_meta[colnames(counts.keep), c("Stage", "sampleNames", "Cell.type", "cross", "Age", "Stage.Age")]
```
make dgelist
```
dgeObj <- DGEList(counts=counts.keep, group=paste0(full_LD_meta$Stage, "_", full_LD_meta$Age))
#add lengths
dgeObj$genes   <- data.frame(Length=edgeR_gene_lengths$Length)
# Can see library size for each sample in the samples slot
dgeObj$samples
```
## Quality control
A few plots to check that data is good quality and that samples are as expected.

Library sizes. Check how many reads per sample in the `dgeObj`.
```dgeObj$samples[,"lib.size"]
```
Distribution. Plot library sizes as barplot to see whether  any major discrepanies between the samples.
```barplot(dgeObj$samples$lib.size, names=colnames(dgeObj), las=2, main="Barplot of library sizes")
abline(h=median(dgeObj$samples$lib.size), lty=2)
```
count data not normally distributed so use log2 scale in rpkm function - corrected for different library sizes.
Can plot boxplots to check distribution of read counts on log2 scale. Add horizontal line corresponding to median logCPM
```rpkm_logcounts <- rpkm(dgeObj,log=TRUE)

boxplot(rpkm_logcounts, xlab="", ylab="Log2 FPKM",las=2)
abline(h=median(rpkm_logcounts), col="blue", main="Boxplots of logFPKMs (unnormalised)")
```
Can also visualise a PCA,which determines the greatest sources of variation in the data, using a multidimensional scaling plot (MDS plot)
colour schemesfor stage (ignoring age for now) 
```full_LD_meta$Stage <- 
  factor(full_LD_meta$Stage)
levels(full_LD_meta$Stage)
col.stage <- c("purple","pink","blue", "lightblue")[full_LD_meta$Stage]
col.stage
```
PlotMDS with stage coloring. Add legend to plot. 
```plotMDS(dgeObj,col=col.stage, main="MDS plot colored by stage", pch=19)
legend("topright",fill=c("purple","pink","blue", "lightblue"),
       legend=levels(full_LD_meta$Stage), bg = "grey", cex = 0.75)
```
## Normalise for composition bias using trimmed mean of M-values normalization method (TMM) by scaling relative to one sample.
Apply normalisation to DGEList object and take a look at normalisation factors for samples
```dgeObj <- calcNormFactors(dgeObj, method="TMM")
dgeObj$samples
head(dgeObj$counts)
```
create an expression matrix to work on if not interested in differential expression - i.e. plotting things.
```ase_fpkm_counts_norm <- rpkm(dgeObj,log=FALSE)
dim(ase_fpkm_counts_norm)
```
can save this if want to use for plotting etc
```saveRDS(ase_fpkm_counts_norm, file = "ase_LD_expression_fpkm.rds")
```
make rownames first column
```
ase_fpkm_counts_norm <- as.data.frame(ase_fpkm_counts_norm)

dim(ase_fpkm_counts_norm)
d <- ase_fpkm_counts_norm
ensembl_gene_id <- rownames(d)
rownames(d) <- NULL
ase_fpkm_counts_norm <- cbind(ensembl_gene_id,d)
dim(ase_fpkm_counts_norm)
```
Melt and combine with metadata and gene annotation to make a data frae for plotting
```
LD_ase_m <- melt(ase_fpkm_counts_norm, id.vars="ensembl_gene_id")
names(LD_ase_m)[2] <- "sampleNames"
names(LD_ase_m)[1] <- "ensembl_gene_id"
```
Annotate with sample details - will take a min to run
```
LD_ase_m <- as.data.frame(LD_ase_m)
annotated <- inner_join(LD_ase_m, full_LD_meta, by="sampleNames")
head(annotated)

#now annotate our data with gene details
ensembl      <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id   <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
                                   'chromosome_name', 'start_position', 'gene_biotype'), mart = ensembl)  
head(ensEMBL2id)
```
Annotation with Gene Details - will take a min to run
```
ann_LD_ase_m <- inner_join(annotated, ensEMBL2id, by="ensembl_gene_id")
head(ann_LD_ase_m)

origin <- c(paste0(ann_LD_ase_m$sampleNames, "_", ann_LD_ase_m$cross))
parental_origin <- substr(origin, start=18, stop=22)

LD_ase_plotting <- cbind(origin, ann_LD_ase_m)
LD_ase_plotting <- cbind(parental_origin, LD_ase_plotting)
LD_ase_plotting$parental_origin <- as.character(LD_ase_plotting$parental_origin)

LD_ase_plotting <- LD_ase_plotting %>% 
  mutate(parent = if_else(parental_origin == "g1_BC" | parental_origin == "g2_CB", "Maternal", "Paternal"))

levels(as.factor(LD_ase_plotting$parent))
```
make subset of data for only the genes you want 
```
ase_plots_data1 <- subset(LD_ase_plotting,  
                          external_gene_name =="Gm10425")

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2 <- data_summary(ase_plots_data1, varname="value", 
                    groupnames=c("external_gene_name", "Stage.Age", "parent"))
```
## Plotting
```
ggplot(df2,
       aes(x=factor(Stage.Age, level = c(
         "Nul.d0","Ges.d5.5","Ges.d9.5","Ges.d14.5","Lac.d5","Lac.d10","Lac.d15","Inv.d1","Inv.d6","Inv.d14")),
         y=value, group=parent, color=parent,
         shape=parent)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~external_gene_name) +
  labs(x = "Stage and age", y = "Expression (FPKM)", color="Parental origin", shape = "Parental origin") +
  ggtitle("Parent-of-origin specific expression") +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2) +
  scale_color_brewer(palette="Set1") +
  theme(
    axis.text.x = element_text(angle=45, hjust=0.5, vjust=0.5),
    title=element_text(size=12, face="bold"),
    axis.title.x = element_text(size=12, face="plain"),
    axis.title.y = element_text(size=12, face="plain"),
    legend.position = "right") +
  annotate("rect", xmin = 1.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           alpha = .3) +
  annotate("rect", xmin = 7.5, xmax = 10.5, ymin = -Inf, ymax = Inf,
           alpha = .3) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           alpha = .1) +
  annotate("rect", xmin = 4.5, xmax = 7.5, ymin = -Inf, ymax = Inf,
           alpha = .1)
           ```
