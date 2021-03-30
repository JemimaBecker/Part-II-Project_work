Some stuff 

have a look at what enrichr terms are significant


using the more stringent adjusted P value <0.05

|Term	|	Overlap	|	P.value	|	Adjusted.P.value	|	Odds.Ratio	|	Combined.Score	|	Genes	|	Pathway	|	UP or DOWN	|
|-------------	|	-------------	|	-------------	|	-------------	|	-------------	|	-------------	|	-------------	|	-------------	|	-------------	|
|mammary gland  lact	|	22/104	|	4.03665E-12	|	3.83481E-10	|	7.931325166	|	208.0831288	|	OLAH;OXTR;PANK3;BTN1A1;RGS16;GLYCAM1;IRX3;CITED4;TRF;CSN1S2A;CSN3;CSN2;ANO1;FABP3;LAO1;SLCO2B1;ELF5;THRSP;MUC15;CSN1S1;LALBA;SLC28A3	|	Mouse_Gene_Atlas	|	DOWN	|
|MP:0004047 abnormal milk composition	|	44805	|	1.76161E-08	|	4.40578E-05	|	20.13739574	|	359.5421853	|	LAO1;BTN1A1;CSN3;PLIN2;THRSP;CSN2;CSN1S1;LALBA;JCHAIN	|	MGI_Mammalian_Phenotype_Level_4_2019	|	DOWN	|
|Inflammatory Response	|	18/200	|	0.000151327	|	0.007112373	|	2.890780975	|	25.42750288	|	IL10;PTGIR;CSF3R;RGS16;EBI3;PTAFR;LIF;KCNA3;TACR1;MEFV;PIK3R5;ADGRE1;P2RX4;CCRL2;PDE4B;NLRP3;LCP2;SLC28A2	|	MGI_SigDB_Hallmark_2020	|	UP	|
|Macrophage markers WP2271	|	44473	|	0.000226935	|	0.025416684	|	19.25361236	|	161.5541411	|	CD86;CD83;LYZ2;RAC2	|	WikiPathways_2019_Mouse	|	DOWN	|
|IL-2/STAT5 Signaling	|	16/199	|	0.001178991	|	0.027706279	|	2.547629147	|	17.17890953	|	IL10;CD86;CD83;RGS16;LIF;PTH1R;NDRG1;AGER;P2RX4;SPP1;BCL2;TNFSF11;PLIN2;SLC39A8;GALM;F2RL2	|	MGI_SigDB_Hallmark_2021	|	UP	|
|regulation of angiogenesis (GO:0045765)	|	13/177	|	2.14034E-05	|	0.038782981	|	4.423125306	|	47.55726748	|	SPARC;SPHK1;PTPRM;EMP2;KLF4;ETS1;HSPG2;RUNX1;SFRP2;STIM1;RRAS;ADGRA2;EPHA1	|	GO_Biological_Processes_2018	|	UP	|
|Head and neck cancer	|	44350	|	0.000112943	|	0.040207735	|	54.69359331	|	497.0896481	|	GALR2;LY6D;LOXL4	|	Jensen_Diseases	|	UP	|
|MP:0005591 decreased vasodilation	|	47300	|	3.87078E-05	|	0.048404159	|	9.222972973	|	93.70050153	|	DHFR;SLC4A7;RGS2;EDNRB;PECAM1;IRS2;KCNN4	|	MGI_Mammalian_Phenotype_Level_4_2020	|	DOWN	|![image](https://user-images.githubusercontent.com/67189202/111080128-ec85a400-84f4-11eb-8f46-9c1d9ee50cca.png)


Having a look at the genes that appear in more than one category:

|	Gene	|	Occurrences	| Notes | chromosome |
|-------------	|	-------------	| -------------	| -------------	|
|	BTN1A1	|	2	| not within 1mp of diff exp lncRNA :(| |
|	CD83	|	2	| within 1mb. chr13 |13|
|	CD86	|	2	|not within 1mp of diff exp lncRNA :(||
|	CSN1S1	|	2	| within 1mb chr 5 |5|
|	CSN2	|	2	| within 1mb chr5 |5|
|	CSN3	|	2	| within 1mb chr5 |5|
|	IL10	|	2	| not within 1mb of lncRNA ||
|	LALBA	|	2	| within 1mb chr15 |15|
|	LAO1	|	2	| within 1mb chr4|4|
|	LIF	|	2	| not within 1mb ||
|	P2RX4	|	2	| within 1mb chr5 |5|
|	PLIN2	|	2	| not within 1mb ||
|	RGS16	|	3	| within 1mb chr1 |1|
|	THRSP	|	2	| not within 1mb ||


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
