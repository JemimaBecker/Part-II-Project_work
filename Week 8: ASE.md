> @jemimabecker also, if you create a column with just sample names (remove the .g1 or .g2) you can filter to the same set as the main expression table. The selection of "clean" samples was done across both datasets. (but works better for the normal expression due to much higher depth). The ASE read counts are MUCH lower, as only a small proportion of reads contain informative SNPs for parental assignment.

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

#return data for genes of interest

```


> Q: is there a good package for looking at this?

> A: "@jemimabecker there is, and I have tried running it, which is how I know the data is still quite noisy
> https://www.bioconductor.org/packages/release/bioc/html/ISoLDE.html
> Actually DESeq2 gives quite good results with a differential across parent genomes"

tutorials and information to work with
- https://bioconductor.org/packages/release/bioc/manuals/ISoLDE/man/ISoLDE.pdf
- https://www.bioconductor.org/packages/release/bioc/html/ISoLDE.html
