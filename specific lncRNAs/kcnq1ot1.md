# Kcnq10t1

import relevant libraries
library("reshape")
library("tidyr")
library("dplyr")
library("ggplot2")

dataframe prep
df <- MouseMammaryGland_Cleaned_edgeR_fpkm_grouped
df.m <- melt(df)
df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("Luminal.Differentiated", "Luminal_differentiated", df.m)))
df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("Luminal.Progenitors", "Luminal_progenitors", df.m)))

df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("14.5", "14_5", df.m)))
df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("5.5", "5_5", df.m)))
df.m <- as.data.frame(lapply(df.m, function(df.m) gsub("9.5", "9_5", df.m)))

df.m.split <- df.m %>% separate(variable, c("Cell.type", "Stage", "Age"), "\\.")

df.m.split <- as.data.frame(lapply(df.m.split, function(df.m) gsub("14_5", "14.5", df.m)))
df.m.split <- as.data.frame(lapply(df.m.split, function(df.m) gsub("5_5", "5.5", df.m)))
df.m.split <- as.data.frame(lapply(df.m.split, function(df.m) gsub("9_5", "9.5", df.m)))

## 1: Endothelial
plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000101609" & Cell.type=="Endothelial"), 
               aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ Stage) +
  theme(axis.text.y=element_blank()) 
plot + labs(x = "Age", y= "Expression (FPKM)", title="Kcnq10t1 expression",ylab=NULL,subtitle = "Endothelial cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Kcnq10t1%20endothelial.png)


## 2: Adipocyte
plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000101609" & Cell.type=="Adipocytes"), 
               aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ Stage) +
  theme(axis.text.y=element_blank()) 
plot + labs(x = "Age", y= "Expression (FPKM)", title="Kcnq10t1 expression",ylab=NULL,subtitle = "Adipocyte cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Kcnq10t1%20adipocyte.png)

## 3: Basal
plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000101609" & Cell.type=="Basal"), 
               aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ Stage) +
  theme(axis.text.y=element_blank()) 
plot + labs(x = "Age", y= "Expression (FPKM)", title="Kcnq10t1 expression",ylab=NULL,subtitle = "Basal cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Kcnq10t1%20basal.png)

## 4: Luminal differentiated
plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000101609" & Cell.type=="Luminal_differentiated"), 
               aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ Stage) +
  theme(axis.text.y=element_blank()) 
plot + labs(x = "Age", y= "Expression (FPKM)", title="Kcnq10t1 expression",ylab=NULL,subtitle = "Luminal differentiated cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Kcnq10t1%20luminal%20differentiated.png)

## 5: Luminal Progenitor
plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000101609" & Cell.type=="Luminal_progenitors"), 
               aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ Stage) +
  theme(axis.text.y=element_blank()) 
plot + labs(x = "Age", y= "Expression (FPKM)", title="Kcnq10t1 expression",ylab=NULL,subtitle = "Luminal progenitor cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Kcnq10t1%20luminal%20progenitors.png)

## 6: Stromal
plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000101609" & Cell.type=="Stromal"), 
               aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), y=value, colour=Stage, group = 3)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ Stage) +
  theme(axis.text.y=element_blank()) 
plot + labs(x = "Age", y= "Expression (FPKM)", title="Kcnq10t1 expression",ylab=NULL,subtitle = "Stromal cells") 

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Kcnq10t1%20stromal.png)
