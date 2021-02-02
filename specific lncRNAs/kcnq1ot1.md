# Kcnq1ot1 (ENSMUSG00000101609)

load libraries

> library("reshape")
> library("tidyr")
> library("dplyr")
> library("ggplot2")

dataframe prep
> df <- MouseMammaryGland_Cleaned_edgeR_fpkm_grouped
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

## Luminal differentiated
> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000101609" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage) 
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Kcnq1ot1 expression",subtitle = "Luminal differentiated cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/kcnq1ot1%20ld.png)

## Luminal Progenitor
> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000101609" & Cell.type=="Luminal_progenitors"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage) +
>   scale_y_continuous(breaks=c(0,1,2,3,4,5,6,7))
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Kcnq1ot1 expression",subtitle = "Luminal progenitor cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/kcnq1otmlp.png)
