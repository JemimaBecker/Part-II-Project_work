# Splitting and orgnaising data
load libraries
> library("reshape")
> library("tidyr")
> library("dplyr")
> library("ggplot2")

start organising grouped data
> df <- MouseMammaryGland_Cleaned_edgeR_fpkm_grouped
> head(df)
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-01-28%20at%2014.32.50.png)
# Reshape the data frame to make it easier to manipulate for plotting
> df.m <- melt(df)
> head(df.m)

## Split the sample name/group into Cell.type, Stage and Age
> df.m.split <- df.m %>% separate(variable, c("Cell.type", "Stage", "Age"), "\\.")
> head(df.m.split)

## Try a simple plot with points
> ggplot(data=subset(df.m.split, X=="ENSMUSG00000000001" & Cell.type=="Endothelial"), 
       aes(x=factor(Age,levels=c("d0","d1","d5","d6","d9","d10","d14","d15")), y=value, colour=Stage)) +
  geom_point())
![] (https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-01-28%20at%2014.39.44.png)


## SPLIT INTO SEPARATE PLOTS FOR STAGE
> ggplot(data=subset(df.m.split, X=="ENSMUSG00000000001" & Cell.type=="Endothelial"), 
       aes(x=factor(Age,levels=c("d0","d1","d5","d6","d9","d10","d14","d15")), y=value, colour=Stage)) +
  geom_point() +
  facet_wrap(~ Stage)
![] (https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-01-28%20at%2014.39.53.png)
