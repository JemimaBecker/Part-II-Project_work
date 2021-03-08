if(!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("IRanges")
library(IRanges)

CTCF_locations <- mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021 %>% dplyr::filter(mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021$V3 ==
                                                                                                            "CTCF_binding_site")
chr1_lncRNA_startstop <- overlaps_chr2_names
ir1 <- IRanges(start=chr1_lncRNA_startstop$end_position,end=chr1_lncRNA_startstop$start_position)

chr1_CTCF <- CTCF_locations %>% dplyr::filter(CTCF_locations$V1 == "2")
ir2 <- IRanges(start =chr1_CTCF$V4,end = chr1_CTCF$V5)

overlaps_table <- table(!is.na(findOverlaps(ir1, ir2, select="arbitrary")))
overlaps_df <- as.data.frame(!is.na(findOverlaps(ir1, ir2, select="arbitrary")))
overlaps_chr2_names <- cbind(overlaps_df,chr1_lncRNA_startstop)
colnames(overlaps_chr2_names)[1] <- "CTCF site"

enhancer_locations <- mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021 %>% dplyr::filter(mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021$V3 ==
                                                                                                            "enhancer")
chr1_lncRNA_startstop <- overlaps_chr2_names
ir1 <- IRanges(start=chr1_lncRNA_startstop$end_position,end=chr1_lncRNA_startstop$start_position)

chr1_enhancers <- enhancer_locations %>% dplyr::filter(enhancer_locations$V1 == "2")
ir2 <- IRanges(start =chr1_enhancers$V4,end = chr1_enhancers$V5)

overlaps_table <- table(!is.na(findOverlaps(ir1, ir2, select="arbitrary")))
overlaps_df <- as.data.frame(!is.na(findOverlaps(ir1, ir2, select="arbitrary")))
overlaps_chr2_names <- cbind(overlaps_df,chr1_lncRNA_startstop)
colnames(overlaps_chr2_names)[1] <- "Enhancer"


promoter_locations <- mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021 %>% dplyr::filter(mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021$V3 ==
                                                                                                            "promoter")
chr1_lncRNA_startstop <- overlaps_chr2_names
ir1 <- IRanges(start=chr1_lncRNA_startstop$end_position,end=chr1_lncRNA_startstop$start_position)

chr1_promoters <- promoter_locations %>% dplyr::filter(promoter_locations$V1 == "2")
ir2 <- IRanges(start =chr1_promoters$V4,end = chr1_promoters$V5)

overlaps_table <- table(!is.na(findOverlaps(ir1, ir2, select="arbitrary")))
overlaps_df <- as.data.frame(!is.na(findOverlaps(ir1, ir2, select="arbitrary")))
overlaps_chr2_names <- cbind(overlaps_df,chr1_lncRNA_startstop)
colnames(overlaps_chr2_names)[1] <- "Promoter"

promoter_flank_locations <- mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021 %>% dplyr::filter(mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021$V3 ==
                                                                                                            "promoter_flanking_region")
chr1_lncRNA_startstop <- overlaps_chr2_names
ir1 <- IRanges(start=chr1_lncRNA_startstop$end_position,end=chr1_lncRNA_startstop$start_position)

chr1_promoterflanks <- promoter_flank_locations %>% dplyr::filter(promoter_flank_locations$V1 == "2")
ir2 <- IRanges(start =chr1_promoterflanks$V4,end = chr1_promoterflanks$V5)

overlaps_table <- table(!is.na(findOverlaps(ir1, ir2, select="arbitrary")))
overlaps_df <- as.data.frame(!is.na(findOverlaps(ir1, ir2, select="arbitrary")))
overlaps_chr2_names <- cbind(overlaps_df,chr1_lncRNA_startstop)
colnames(overlaps_chr2_names)[1] <- "Promoter flanking region"


TF_locations <- mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021 %>% dplyr::filter(mus_musculus.GRCm39.Regulatory_Build.regulatory_features.20201021$V3 ==
                                                                                                            "TF_binding_site")
chr1_lncRNA_startstop <- overlaps_chr2_names
ir1 <- IRanges(start=chr1_lncRNA_startstop$end_position,end=chr1_lncRNA_startstop$start_position)

chr1_TFs <- TF_locations %>% dplyr::filter(TF_locations$V1 == "2")
ir2 <- IRanges(start =chr1_TFs$V4,end = chr1_TFs$V5)

overlaps_table <- table(!is.na(findOverlaps(ir1, ir2, select="arbitrary")))
overlaps_df <- as.data.frame(!is.na(findOverlaps(ir1, ir2, select="arbitrary")))
overlaps_chr2_names <- cbind(overlaps_df,chr1_lncRNA_startstop)
colnames(overlaps_chr2_names)[1] <- "TF Binding site"

