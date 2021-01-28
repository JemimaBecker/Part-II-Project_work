# Examination of the imprinted miRNA containing gene Mirg
## Initial data organisation
> mirg.dataG <- ncExpG %>% dplyr::filter(ncExpG$external_gene_name == "Mirg")
> mirg.dataG <- t(mirg.dataG)
> colnames(mirg.dataG) = mirg.dataG[1, ] 

### ANOVA: Does the data fit the requirements and assumptions of this model?
> boxplot(exp ~ Cell.type, data=mirg)
> mirg.lm <- lm(exp ~ Cell.type, data=mirg)
> anova(mirg.lm)
> par(mfrow=c(2,2))
> plot(mirg.lm)
data doesnt seem to fit the anova criteria very well - normal QQ plot is particularly bad
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-01-28%20at%2014.51.33.png)
> kruskal.test(exp ~ Cell.type, data=mirg)
significant variance in mirg expression between cell types
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-01-28%20at%2014.51.03.png)
> kruskal.test(exp ~ Dev.state, data=mirg)
not significant variant in exp by developmental state, but will look at individual cell types
> par(mfrow=c(1,1))
> PlotViolin(exp ~ Cell.type, data=mirg)
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-01-28%20at%2014.51.53.png)
# mirg in various cell types:

## 1: stroma
> mirg.stroma <- mirg %>% dplyr::filter(mirg$Cell.type == "Stromal")
> boxplot(exp ~ Dev.state, data=mirg.stroma)
> kruskal.test(exp ~ Dev.state, data=mirg.stroma)
data:  exp by Dev.state
Kruskal-Wallis chi-squared = 0.92727, df = 3, p-value = 0.8188
> mirg.lm <- lm(exp ~ Cell.type, data=mirg)
> anova(mirg.lm)
> par(mfrow=c(2,2))
> plot(mirg.lm)

> par(mfrow=c(2,3))
> mirg.stroma <- mirg %>% dplyr::filter(mirg$Cell.type == "Stromal")
> boxplot(exp ~ Dev.state, data=mirg.stroma)

## 2: Adipocytes
> mirg.adi <- mirg %>% dplyr::filter(mirg$Cell.type == "Adipocyte")
> boxplot(exp ~ Dev.state, data=mirg.adi, main="Mirg exp in adipocytes")
> kruskal.test(exp ~ Dev.state, data=mirg.adi)
data:  exp by Dev.state
Kruskal-Wallis chi-squared = 2.4621, df = 3, p-value = 0.4822

## 3: Basal
> mirg.basal <- mirg %>% dplyr::filter(mirg$Cell.type == "Basal")
> boxplot(exp ~ Dev.state, data=mirg.basal, main="Mirg exp in basal cells")
> kruskal.test(exp ~ Dev.state, data=mirg.basal)
data:  exp by Dev.state
Kruskal-Wallis chi-squared = 4.1707, df = 3, p-value = 0.2436

## 4: Endothelial cells
> mirg.endo <- mirg %>% dplyr::filter(mirg$Cell.type == "Endothelial")
> boxplot(exp ~ Dev.state, data=mirg.endo, main="Mirg exp in endothelial cells")
> kruskal.test(exp ~ Dev.state, data=mirg.endo)
data:  exp by Dev.state
Kruskal-Wallis chi-squared = 1.061, df = 3, p-value = 0.7865

## 5: luminal differentiated cells
> mirg.LD <- mirg %>% dplyr::filter(mirg$Cell.type == "Luminal.differentiated")
> boxplot(exp ~ Dev.state, data=mirg.LD, main="Mirg exp in LD cells")
> kruskal.test(exp ~ Dev.state, data=mirg.LD)
data:  exp by Dev.state
Kruskal-Wallis chi-squared = 5.1455, df = 3, p-value = 0.1615


## 6: luminal progenitor cells
> mirg.LP <- mirg %>% dplyr::filter(mirg$Cell.type == "Luminal.progenitor")
> boxplot(exp ~ Dev.state, data=mirg.LP, main="Mirg exp in LP cells")
> kruskal.test(exp ~ Dev.state, data=mirg.LP)
data:  exp by Dev.state
Kruskal-Wallis chi-squared = 0.40244, df = 3, p-value = 0.9397
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-01-28%20at%2014.52.21.png)
# Mirg at various time points

## 1 Nulliparous
> mirg.NP <- mirg %>% dplyr::filter(mirg$Dev.state == "Nulliparous")
> boxplot(exp ~ Cell.type, data=mirg.NP, main="Mirg exp Nulliparous")
> kruskal.test(exp ~ Cell.type, data=mirg.NP)
data:  exp by Cell.type
Kruskal-Wallis chi-squared = 5, df = 5, p-value = 0.4159

## 2 Gestation
> mirg.Gest <- mirg %>% dplyr::filter(mirg$Dev.state == "Gestation")
> boxplot(exp ~ Cell.type, data=mirg.Gest, main="Mirg exp Gestation")
> kruskal.test(exp ~ Cell.type, data=mirg.Gest)
data:  exp by Cell.type
Kruskal-Wallis chi-squared = 6.334, df = 5, p-value = 0.2751

## 3 Lactation
> mirg.Lact <- mirg %>% dplyr::filter(mirg$Dev.state == "Lactation")
> boxplot(exp ~ Cell.type, data=mirg.Lact, main="Mirg exp Lactation")
> kruskal.test(exp ~ Cell.type, data=mirg.Lact)
data:  exp by Cell.type
Kruskal-Wallis chi-squared = 6.7007, df = 5, p-value = 0.2439

## 4 Involution
> mirg.INV <- mirg %>% dplyr::filter(mirg$Dev.state == "Involution")
> boxplot(exp ~ Cell.type, data=mirg.INV, main="Mirg exp Involution")
> kruskal.test(exp ~ Cell.type, data=mirg.INV)
data:  exp by Cell.type
Kruskal-Wallis chi-squared = 7.6385, df = 5, p-value = 0.1773
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-01-28%20at%2014.52.55.png)
