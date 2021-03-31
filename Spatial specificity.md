```
imp_lnc_plot <- df.m.split %>% dplyr::filter(df.m.split$X=="ENSMUSG00000078247"|
X=="ENSMUSG00000086537"|
X=="ENSMUSG00000101609"|
X=="ENSMUSG00000021268" |
X=="ENSMUSG00000100826" |
X=="ENSMUSG00000000031")

ggplot(data=subset(imp_lnc_plot, X=="ENSMUSG00000100826" & Stage=="Gestation"), 
       aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
           y=value, colour=Cell.type, group = paste0(Stage, "_", Cell.type))) +
  geom_point() +
  geom_line() +
  facet_wrap(~ X)  +
  labs(x = "Days (from start of stage)", y= "Expression (FPKM)", title="Gestation expression",subtitle = "All Cell Types") +
  theme_bw() +
  theme(legend.position="right")
```



```
red_nulliparous <- read.csv("/data/homes/jb2220/red_nulliparous.csv")
red_gestation <- read.csv("/data/homes/jb2220/red_gestation.csv")
red_lactation <- read.csv("/data/homes/jb2220/red_lactation.csv")
red_involution <- read.csv("/data/homes/jb2220/red_involution.csv")

rownames(red_nulliparous) <- red_nulliparous[,1]
rownames(red_gestation) <- red_gestation[,1]
rownames(red_lactation) <- red_lactation[,1]
rownames(red_involution) <- red_involution[,1]


ensembl    <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description',
                                 'chromosome_name', 'start_position', 'end_position', 
                                 'strand', 'gene_biotype'), mart = ensembl)  
rownames(ensEMBL2id) <- ensEMBL2id$ensembl_gene_id

red_ann_nulliparous  <- merge(red_nulliparous, ensEMBL2id, by=0)
red_ann_gestation  <- merge(red_gestation, ensEMBL2id, by=0)
red_ann_lactation  <- merge(red_lactation, ensEMBL2id, by=0)
red_ann_involution  <- merge(red_involution, ensEMBL2id, by=0)


red_ann_NC_nulliparous <- red_ann_nulliparous %>% dplyr::filter(red_ann_nulliparous$gene_biotype == 
                                                   "3prime_overlapping_ncRNA" |
                                                   gene_biotype == "antisense" |
                                                   gene_biotype == "bidirectional_promoter_lncRNA" |
                                                   gene_biotype == "lincRNA" |
                                                   gene_biotype == "lncRNA" |
                                                   gene_biotype == "macro_lncRNA" |
                                                   gene_biotype ==  "processed_transcript")

red_ann_NC_gestation <- red_ann_gestation %>% dplyr::filter(red_ann_gestation$gene_biotype == 
                                                             "3prime_overlapping_ncRNA" |
                                                             gene_biotype == "antisense" |
                                                             gene_biotype == "bidirectional_promoter_lncRNA" |
                                                             gene_biotype == "lincRNA" |
                                                             gene_biotype == "lncRNA" |
                                                             gene_biotype == "macro_lncRNA" |
                                                             gene_biotype ==  "processed_transcript")
red_ann_NC_lactation <- red_ann_lactation %>% dplyr::filter(red_ann_lactation$gene_biotype == 
                                                             "3prime_overlapping_ncRNA" |
                                                             gene_biotype == "antisense" |
                                                             gene_biotype == "bidirectional_promoter_lncRNA" |
                                                             gene_biotype == "lincRNA" |
                                                             gene_biotype == "lncRNA" |
                                                             gene_biotype == "macro_lncRNA" |
                                                             gene_biotype ==  "processed_transcript")
red_ann_NC_involution<- red_ann_involution %>% dplyr::filter(red_ann_involution$gene_biotype == 
                                                             "3prime_overlapping_ncRNA" |
                                                             gene_biotype == "antisense" |
                                                             gene_biotype == "bidirectional_promoter_lncRNA" |
                                                             gene_biotype == "lincRNA" |
                                                             gene_biotype == "lncRNA" |
                                                             gene_biotype == "macro_lncRNA" |
                                                             gene_biotype ==  "processed_transcript")


spatial_master <- rbind(red_ann_NC_nulliparous,red_ann_NC_gestation,red_ann_NC_lactation,red_ann_NC_involution)
```
