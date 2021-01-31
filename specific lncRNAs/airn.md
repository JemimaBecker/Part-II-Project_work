current state (31-01-21):
need to find a way to split data without splitting luminal.progenitor into two columns. current work bypasses luminal samples until code can be modified

## load packages
imprinted_lncRNAs <- c("ENSMUSG00000078247","ENSMUSG00000101609","ENSMUSG00000086537","ENSMUSG00000100826","ENSMUSG00000000031","ENSMUSG00000021268")

library("reshape")
library("tidyr")
library("dplyr")
library("ggplot2")

## dataframe prep
df <- MouseMammaryGland_Cleaned_edgeR_fpkm_grouped
df.m <- melt(df)
df.m.split <- df.m %>% separate(variable, c("Cell.type", "Stage", "Age"), "\\.")
## needs to be modified so that Luminal.Differentiated and Luminal.Progenitors arent split into two separate columns


# airn
## plot by time in four different stages for specific cell type
1: Endothelial
> ggplot(data=subset(df.m.split, X=="ENSMUSG00000078247" & Cell.type=="Endothelial"), 
>        aes(x=factor(Age,levels=c("d0","d1","d5","d6","d9","d10","d14","d15")), y=value, colour=Stage, group = 2)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage)
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-01-31%20at%2016.53.18.png)
2: Adipocyte

> ggplot(data=subset(df.m.split, X=="ENSMUSG00000078247" & Cell.type=="Adipocytes"), 
>        aes(x=factor(Age,levels=c("d0","d1","d5","d6","d9","d10","d14","d15")), y=value, colour=Stage, group = 2)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage)

3: Basal

> ggplot(data=subset(df.m.split, X=="ENSMUSG00000078247" & Cell.type=="Basal"), 
>        aes(x=factor(Age,levels=c("d0","d1","d5","d6","d9","d10","d14","d15")), y=value, colour=Stage, group = 2)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage)

4: Luminal differentiated

5: Luminal Progenitor

6: Stromal

> ggplot(data=subset(df.m.split, X=="ENSMUSG00000078247" & Cell.type=="Stromal"), 
>        aes(x=factor(Age,levels=c("d0","d1","d5","d6","d9","d10","d14","d15")), y=value, colour=Stage, group = 2)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage)
