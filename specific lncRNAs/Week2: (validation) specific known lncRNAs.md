> imprinted_lncRNAs <- c("ENSMUSG00000078247","ENSMUSG00000101609","ENSMUSG00000086537","ENSMUSG00000100826","ENSMUSG00000000031","ENSMUSG00000021268")

> imprinted_known_FPKMG <- MouseMammaryGland_Cleaned_edgeR_fpkm_grouped %>% 
>   dplyr::filter(MouseMammaryGland_Cleaned_edgeR_fpkm_grouped$X == "ENSMUSG00000078247" | 
>            X == "ENSMUSG00000101609" | 
>            X == "ENSMUSG00000086537" | 
>            X == "ENSMUSG00000100826" | 
>            X == "ENSMUSG00000000031" | 
>            X == "ENSMUSG00000021268")
> library("reshape")
> library("tidyr")
> library("dplyr")
> library("ggplot2")
## dataframe prep
> df <- imprinted_known_FPKMG
> df.m <- melt(df)
> df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("Luminal.Differentiated", "Luminal_differentiated", df.m)))
> df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("Luminal.Progenitors", "Luminal_progenitors", df.m)))

> df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("d14.5", "d14_5", df.m)))
> df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("d5.5", "d5_5", df.m)))
> df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("d9.5", "d9_5", df.m)))

> df.m.split <- df.m %>% separate(variable, c("Cell.type", "Stage", "Age"), "\\.")

> df.m.split <- as.data.frame(lapply(df.m.split, function(df.m) gsub("d14_5", "d14.5", df.m)))
> df.m.split <- as.data.frame(lapply(df.m.split, function(df.m) gsub("d5_5", "d5.5", df.m)))
> df.m.split <- as.data.frame(lapply(df.m.split, function(df.m) gsub("d9_5", "d9.5", df.m)))

> df.m.split$value <- as.numeric(as.character(df.m.split$value))

## Work to be done

need to look now at significance of variation of expression

- [ ] Q1: is developmental stage a strong predictor of expression? (for each individual)

- [ ]  Q2: is there signifiant varition within each stage (of expression per individual)

> i.e is age (d) a significant predictor of expression within each stage?

- [ ] Q3: is there coordinate regulation of genes within 1Mb in cis from any of these that show significant differences in stage-specific expression?
