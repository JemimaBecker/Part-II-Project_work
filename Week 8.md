Some stuff 

have a look at what enrichr terms are significant

```
library(devtools)
install_github("wjawaid/enrichR")
install.packages("enrichR")
library(enrichR)
listEnrichrSites()
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)
dbs <- c(
  "Allen_Brain_Atlas_10x_scRNA_2021",
  "Allen_Brain_Atlas_down",
  "Allen_Brain_Atlas_up",
  "BioPlanet_2019",
  "GO_Biological_Process_2018",
  "GO_Cellular_Component_2018",
  "GO_Molecular_Function_2018",
  "Elsevier_Pathway_Collection",
  "MGI_Mammalian_Phenotype_Level_4_2019",
  "MSigDB_Hallmark_2020",
  "Mouse_Gene_Atlas",
  "Jensen_DISEASES",
  "KEGG_2019_Mouse",
  "WikiPathways_2019_Mouse",
  "Reactome_2016")
if (websiteLive) {
  enriched <- enrichr(c(**LIST OF DOWNREGULATED GENES***), dbs)
}

UP_Allen_Brain_Atlas_down <- as.data.frame(if (websiteLive) enriched[["Allen_Brain_Atlas_down"]])
UP_Allen_Brain_Atlas_up <- as.data.frame(if (websiteLive) enriched[["Allen_Brain_Atlas_up"]])
UP_BioPlanet_2019 <- as.data.frame(if (websiteLive) enriched[["BioPlanet_2019"]])
UP_GO_Biological_Process_2018 <- as.data.frame(if (websiteLive) enriched[["GO_Biological_Process_2018"]])
UP_GO_Cellular_Component_2018 <- as.data.frame(if (websiteLive) enriched[["GO_Cellular_Component_2018"]])
UP_GO_Molecular_Function_2018 <- as.data.frame(if (websiteLive) enriched[["GO_Molecular_Function_2018"]])
UP_Elsevier_Pathway_Collection <- as.data.frame(if (websiteLive) enriched[["Elsevier_Pathway_Collection"]])
UP_MGI_Mammalian_Phenotype_Level_4_2019 <- as.data.frame(if (websiteLive) enriched[["MGI_Mammalian_Phenotype_Level_4_2019"]])
UP_MSigDB_Hallmark_2020 <- as.data.frame(if (websiteLive) enriched[["MSigDB_Hallmark_2020"]])
UP_Mouse_Gene_Atlas <- as.data.frame(if (websiteLive) enriched[["Mouse_Gene_Atlas"]])
UP_Jensen_DISEASES <- as.data.frame(if (websiteLive) enriched[["Jensen_DISEASES"]])
UP_KEGG_2019_Mouse <- as.data.frame(if (websiteLive) enriched[["KEGG_2019_Mouse"]])
UP_WikiPathways_2019_Mouse <- as.data.frame(if (websiteLive) enriched[["WikiPathways_2019_Mouse"]])
UP_Reactome_2016 <- as.data.frame(if (websiteLive) enriched[["Reactome_2016"]])

write.csv(UP_Allen_Brain_Atlas_down,file="UP_Allen_Brain_Atlas_down.csv")
write.csv(UP_Allen_Brain_Atlas_up,file="UP_Allen_Brain_Atlas_up.csv")
write.csv(UP_BioPlanet_2019,file="UP_BioPlanet_2019.csv")
write.csv(UP_GO_Biological_Process_2018,file="UP_GO_Biological_Process_2018.csv")
write.csv(UP_GO_Cellular_Component_2018,file="UP_GO_Cellular_Component_2018.csv")
write.csv(UP_GO_Molecular_Function_2018,file="UP_GO_Molecular_Function_2018.csv")
write.csv(UP_Elsevier_Pathway_Collection,file="UP_Elsevier_Pathway_Collection.csv")
write.csv(UP_MGI_Mammalian_Phenotype_Level_4_2019,file="UP_MGI_Mammalian_Phenotype_Level_4_2019.csv")
write.csv(UP_MSigDB_Hallmark_2020,file="UP_MSigDB_Hallmark_2020.csv")
write.csv(UP_Mouse_Gene_Atlas,file="UP_Mouse_Gene_Atlas.csv")
write.csv(UP_Jensen_DISEASES,file="UP_Jensen_DISEASES.csv")
write.csv(UP_KEGG_2019_Mouse,file="UP_KEGG_2019_Mouse.csv")
write.csv(UP_WikiPathways_2019_Mouse,file="UP_WikiPathways_2019_Mouse.csv")
write.csv(UP_Reactome_2016,file="UP_Reactome_2016.csv")
```
