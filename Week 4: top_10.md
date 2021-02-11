# Analysis of top 10 differentially expressed lncRNAs in differentiated luminal epithelium cells

These are the 10 genes with the lowest adjusted p-values as returned by DESeq2 (top 10 from annotated_red.csv)

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-11%20at%2018.22.37.png)

## 1: Prepare data

load libraries

> library("reshape")
> library("tidyr")
> library("dplyr")
> library("ggplot2")

Reformat FPKM data

> par(mfrow=c(2,2))
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
> df.m.split$Stage_ordered <- factor(df.m.split$Stage, 
>                              levels=c('Nulliparous','Gestation','Lactation','Involution'))
> saveRDS(df.m.split, file = "df.m.split.rds")

## 2: Plot gene expression by stage

### 1:ENSMUSG00000091423 Gm17509

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000091423" & Cell.type=="Luminal_differentiated"), 
>                 aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                     y=value, colour=Stage_ordered, group = 3)) +
>    geom_point() +
>    geom_line() +
>    facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm17509 expression (ENSMUSG00000091423)",subtitle = "Luminal differentiated cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm17509.png)

### 2:ENSMUSG00000085649 A730032A03Rik

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000085649" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="A730032A03Rik expression (ENSMUSG00000085649)",subtitle = "Luminal differentiated cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/a730032a03rik.png)

### 3:ENSMUSG00000089961 Gm16567

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000089961" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm16567 expression (ENSMUSG00000089961)",subtitle = "Luminal differentiated cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm16567.png)

### 4:ENSMUSG00000100954 Gm10138

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000100954" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm10138 expression (ENSMUSG00000100954)",subtitle = "Luminal differentiated cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm10138.png)

### 5:ENSMUSG00000090208 Gm15851

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000090208" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm15851 expression (ENSMUSG00000090208)",subtitle = "Luminal differentiated cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm15851.png)

### 6:ENSMUSG00000087698 Gm13031

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000087698" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm13031 expression (ENSMUSG00000087698)",subtitle = "Luminal differentiated cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm13031.png)

### 7:ENSMUSG00000086052 Gm11802

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000086052" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm11802 expression (ENSMUSG00000086052)",subtitle = "Luminal differentiated cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm11802.png)

### 8:ENSMUSG00000092283 Gm20412

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000092283" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm20412 expression (ENSMUSG00000092283)",subtitle = "Luminal differentiated cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm20412.png)

### 9:ENSMUSG00000097988 Gm10535

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000097988" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm10535 expression (ENSMUSG00000097988)",subtitle = "Luminal differentiated cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm10535.png)

### 10:ENSMUSG00000085083 Gm11615

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000085083" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm11615 expression (ENSMUSG00000085083)",subtitle = "Luminal differentiated cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm11615.png)
