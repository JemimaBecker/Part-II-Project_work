# Airn (ENSMUSG00000078247)
current state (31-01-21):
need to find a way to split data without splitting luminal.progenitor into two columns. current work bypasses luminal samples until code can be modified
update (01-02-210
should be sorted now

## load packages
> library("reshape")
> library("tidyr")
> library("dplyr")
> library("ggplot2")

## dataframe prep
> df <- MouseMammaryGland_Cleaned_edgeR_fpkm_grouped
> df.m <- melt(df)
> df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("Luminal.Differentiated", "Luminal_differentiated", df.m)))
> df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("Luminal.Progenitors", "Luminal_progenitors", df.m)))

> df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("14.5", "14_5", df.m)))
> df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("5.5", "5_5", df.m)))
> df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("9.5", "9_5", df.m)))

> df.m.split <- df.m %>% separate(variable, c("Cell.type", "Stage", "Age"), "\\.")

> df.m.split <- as.data.frame(lapply(df.m.split, function(df.m) gsub("14_5", "14.5", df.m)))
> df.m.split <- as.data.frame(lapply(df.m.split, function(df.m) gsub("5_5", "5.5", df.m)))
> df.m.split <- as.data.frame(lapply(df.m.split, function(df.m) gsub("9_5", "9.5", df.m)))

# airn
## plot by time in four different stages for specific cell type

1: Endothelial

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000078247" & Cell.type=="Endothelial"), 
>        aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage) +
>   theme(axis.text.y=element_blank()) 
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Airn expression",ylab=NULL,subtitle = "Endothelial cells") 

2: Adipocyte

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000078247" & Cell.type=="Adipocytes"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage) +
>   theme(axis.text.y=element_blank()) 
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Airn expression",ylab=NULL,subtitle = "Adipoctye cells") 

3: Basal

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000078247" & Cell.type=="Basal"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage) +
>   theme(axis.text.y=element_blank()) 
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Airn expression",ylab=NULL,subtitle = "Basal cells") 

4: Luminal differentiated

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000078247" & Cell.type=="Luminal_differentiated"), 
>              aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage) +
> theme(axis.text.y=element_blank()) 
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Airn expression",ylab=NULL,subtitle = "Luminal differentiated cells cells") 

5: Luminal Progenitor

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000078247" & Cell.type=="Luminal_progenitors"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage) +
>   theme(axis.text.y=element_blank()) 
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Airn expression",ylab=NULL,subtitle = "Luminal progenitor cells") 

6: Stromal

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000078247" & Cell.type=="Stromal"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage) +
>   theme(axis.text.y=element_blank()) 
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Airn expression",ylab=NULL,subtitle = "Stromal cells") 

