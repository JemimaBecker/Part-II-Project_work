# Analysis of top 10 differentially expressed lncRNAs in differentiated luminal epithelium cells

These are the 10 genes with the lowest adjusted p-values as returned by DESeq2 (top 10 from annotated_red.csv)

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-11%20at%2018.22.37.png)

of these, 9/10 overlap with a CTCF binding site, and all 10 overlap a promoter flank or a promoter.

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

This gene is located on chromosome 13 and is antisense to Embigin (ENSMUSG00000021728), a gene with known roles in developmental control and cell differentiation.
Embigin shows very similar patterns of expression, also its expression varies over a smaller order of magnitude. Futher investigation is needed to examine the causal/functional relationship here.

Embigin padj: 0.001311741, log2FoldChange: -3.861413 (red)

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-12%20at%2017.14.09.png)

> http://www.informatics.jax.org/gxd/marker/MGI:95321?tab=stagegridtab#gxd=markerMgiId%3DMGI%3A95321%26theilerStage%3D%26assayType%3D%26results%3D100%26startIndex%3D0%26sort%3D%26dir%3Dasc%26tab%3Dstagegridtab
> https://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=ENSMUSG00000021728;r=13:117218701-117221075;t=ENSMUST00000022242

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000091423" & Cell.type=="Luminal_differentiated"), 
>                 aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                     y=value, colour=Stage_ordered, group = 3)) +
>    geom_point() +
>    geom_line() +
>    facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm17509 expression (ENSMUSG00000091423)",subtitle = "Luminal differentiated cells") 


|   |  |
| ------------- | ------------- |
| ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm17509.png)|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-12%20at%2017.27.43.png)|



### 2:ENSMUSG00000085649 A730032A03Rik

Antisense to Wfdc5 (ENSMUSG00000040154) - an extracellular protease inhibitor. Weird patterns, doesnt make it through to ither red or bue significance categories

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-12%20at%2017.29.44.png)

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000085649" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="A730032A03Rik expression (ENSMUSG00000085649)",subtitle = "Luminal differentiated cells") 

|   |  |
| ------------- | ------------- |
|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/a730032a03rik.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Screenshot%202021-02-12%20at%2017.36.26.png)|

### 3:ENSMUSG00000089961 Gm16567

Antisense to C1s1 (ENSMUSG00000038521) doesnt make it through to either significance category

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-12%20at%2017.39.41.png)

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000089961" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm16567 expression (ENSMUSG00000089961)",subtitle = "Luminal differentiated cells") 

|   |  |
| ------------- | ------------- |
|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm16567.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-12%20at%2017.41.58.png)|

### 4:ENSMUSG00000100954 Gm10138

Antisense to Ivns1abp (ENSMUSG00000023150) - role in cell death?

Ivns1abp padj: 0.002807668, log2FoldChange: -1.770356 (blue)

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-12%20at%2017.49.46.png)

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000100954" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm10138 expression (ENSMUSG00000100954)",subtitle = "Luminal differentiated cells") 

|   |  |
| ------------- | ------------- |
|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm10138.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-12%20at%2017.51.46.png)

### 5:ENSMUSG00000090208 Gm15851

Antisense to two protein coding genes Optc (ENSMUSG00000010311) and Prelp (ENSMUSG00000041577). Neither make it through to significance categories

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm15851%20location.png)

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000090208" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm15851 expression (ENSMUSG00000090208)",subtitle = "Luminal differentiated cells") 

|   |  |   |
| ------------- | ------------- | ------------- |
| ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm15851.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/optc.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/prelp.png) |

### 6:ENSMUSG00000087698 Gm13031

Antisense to Padi2 (ENSMUSG00000028927, Protein-arginine deiminase type-2) and very close to Sdhb (ENSMUSG00000009863, succinate dehydrogenase complex subunit B)

Padi2: padj: 3.567672e-10, log2FoldChange: 3.935478 (red)

Sdhb variance isnt significant.

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2012.12.50.png)

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000087698" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm13031 expression (ENSMUSG00000087698)",subtitle = "Luminal differentiated cells") 

|   |  |   |
| ------------- | ------------- | ------------- |
|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm13031.png)| ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2012.17.38.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2012.18.15.png) |

### 7:ENSMUSG00000086052 Gm11802

Antisense to Tox (ENSMUSG00000041272, thymocyte selection-associated high mobility group box).  Not significant variance :(

"TOX is also a member of a small subfamily of proteins (TOX2, TOX3, and TOX4) that share almost identical HMG-box sequences.[8] TOX3 has been identified as a breast cancer susceptibility locus." https://en.wikipedia.org/wiki/TOX

Development of all CD4 T lineages requires nuclear factor TOX https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2234360/

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2012.22.18.png)

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000086052" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm11802 expression (ENSMUSG00000086052)",subtitle = "Luminal differentiated cells") 

|   |  |
| ------------- | ------------- |
|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm11802.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2012.26.14.png) |

### 8:ENSMUSG00000092283 Gm20412

within the intron of a Rprd1b (ENSMUSG00000027651)(as is another lncRNA: 2010009K17Rik/ENSMUSG00000100860, not significant change) and antisense to Tgm2 (ENSMUSG00000037820)

Rprd1b promotes cell proliferation and is upregulated in tumours. Padj= 1.914443e-11, log2FoldChange= 3.037133 (red)

upregulation of Tgm2 associated with cancer metastasis and lower survival rates, implicated also in breast cancer. Padj=5.14246 e-05, log2FoldChange= 2.990424 (red)

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2012.32.27.png)

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000092283" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm20412 expression (ENSMUSG00000092283)",subtitle = "Luminal differentiated cells") 

|   |  |
| ------------- | ------------- |
| ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm20412.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2012.43.48.png) |
| ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2012.44.31.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2012.45.05.png) |

### 9:ENSMUSG00000097988 Gm10535

Antisense to Ptprv (ENSMUSG00000097993, not significant) and nearby to Lgr6 (ENSMUSG00000042793, not significant)

Ptprv May play a role in the maintenance of pluripotency. Down-regulated during differentiation. https://www.uniprot.org/uniprot/P70289

Lgr6 "Receptor for R-spondins that potentiates the canonical Wnt signaling pathway and acts as a marker of multipotent stem cells in the epidermis..May act as a tumor suppressor." https://www.uniprot.org/uniprot/Q9HBX8

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2012.53.08.png)

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000097988" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm10535 expression (ENSMUSG00000097988)",subtitle = "Luminal differentiated cells") 

|   |  | |
| ------------- | ------------- | ------------- |
|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm10535.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2012.56.58.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2012.57.37.png) |

### 10:ENSMUSG00000085083 Gm11615

lots of protein coding genes in the near vicinity
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2013.02.31.png) 
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2013.03.53.png)

- Wnk4 ENSMUSG00000035112, padj=1.527823e-06, log2FoldChange= -2.512774 (red). Serine/threonine protein kinase,WNK4 plays a critical role in the regulation of various transporters and channels in the kidney (https://en.wikipedia.org/wiki/WNK4)
- Coa3 ENSMUSG00000017188, padj= 0.02983283, log2FoldChange= -2.821824 (red). Cytochrome c oxidase assembly factor 3
- Cntd1 ENSMUSG00000078653, padj= 0.005747596, log2FoldChange= -0.6838967 (bue). Cyclin N-terminal domain-containing protein 1
- Vps25 ENSMUSG00000078656, not significant. Vacuolar protein-sorting-associated protein 25, Component of the ESCRT-II complex
- Becn1 ENSMUSG00000035086, padj= 8.836736e-06, log2FoldChange= -2.426348 (red). role in regulation of tumourigenesis and cell death
- Ramp2 ENSMUSG00000001240, not significant. Receptor activity modifying protein 2, role in glycosylation and transportation of adrenomedullin receptor to the cell surface

> plot <- ggplot(data=subset(df.m.split, X=="ENSMUSG00000085083" & Cell.type=="Luminal_differentiated"), 
>                aes(x=factor(Age,levels=c("d0","d1","d5","d5.5","d6","d9","d9.5","d10","d14","d14.5","d15")), 
>                    y=value, colour=Stage_ordered, group = 3)) +
>   geom_point() +
>   geom_line() +
>   facet_wrap(~ Stage_ordered, ncol=2 )  
> plot + labs(x = "Age", y= "Expression (FPKM)", title="Gm11615 expression (ENSMUSG00000085083)",subtitle = "Luminal differentiated cells") 

|   |  |
| ------------- | ------------- |
|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/gm11615.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2014.33.29.png) |



|   |  | |
| ------------- | ------------- | ------------- |
|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2014.20.45.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2014.22.04.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2014.21.46.png) |
| ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2014.21.30.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-13%20at%2014.20.55.png) |  |
