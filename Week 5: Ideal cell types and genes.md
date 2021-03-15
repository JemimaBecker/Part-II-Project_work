# Week 5: Identifying the best cell type and candidate genes

## DESeq2

I ran DESeq2 on the lncRNA genes in each cell type individually to examine:

- How many genes were differentially expressed
- What the distribution of expression looked like

| Cell type  |Number of lncRNA genes with padj<0.1 and log2FoldChange > 2| total # diff exp genes |
| ------------- | ------------- |------------- |
|Luminal Differentiated| 234 |840|
|Stromal|132|931|
|Adipocyte|66|217|
|Luminal Progenitor|62|144|
|Basal|40|162|
|Endothelial|23|172|

although the absolute numbers of differentially expressed genes in stromal is higher than in LD, based on the distributions below, LD will be more itneresting.

|||
| ------------- | ------------- |
|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/week%208%20images/Screenshot%202021-03-15%20at%2011.08.40.pnghttps://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/week%208%20images/Screenshot%202021-03-15%20at%2011.11.49.png


Looking at distributions of lncRNA DE in differnt cell types
|||
| ------------- | ------------- |
| ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/luminal_differentiated_nc.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/stromal_nc.png) |
| ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/adipocytes_nc.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/luminal_progenitor_nc.png) |
| ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/basal_nc.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/endothelial_nc.png) |

By comparing the lists of differentially expressed lncRNAs (padj < 0.1, |log2FoldChange|>2) across cell types I obtained a list of 77 lncRNA genes that were significantly differentially expressed in two or more cell types

A small number(11) were differentially expressed in three or more cell types

| Gene name  |Number of cell types | Luminal Differentiated | Luminal Progenitors | Stromal | Basal | Adipoctye | Endothelial |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| Gm15527	| 4| +| + | + | + | -|- |	
| Gm15629	| 4| + |- | + | + | +|- |	
| Gm42697	| 4|  +| + |-  | + | -| + |	
| Gm13031	| 3| + | + | + |- |- |- |	
| Gm14091	| 3| + |- | + |- | -|  +|	
| Gm15650	| 3| + |- |- | + |- |- |	
| Gm28513	| 3| + | + |- | -|- |- |	
| Gm42705	| 3| + | + |- | + | -|- |	
| 8030451A03Rik	| 3| + | + |- | + | |- |	
| Gm10570	| 3| + | + |- | -| -|- |	
| Gimap1os	| 3| + | + |- | -|- | -|	

| Gene name | Differnetially expressed protein coding genes within 1mb |
| ------------- | ------------- |
| Gm15527	| none |
| Gm15629	| Rubcnl, Lcp1|
| Gm42697	| Pla2g12a|
| Gm13031	| Padi2, Epha2,Fblim2,Efhd2|
| Gm14091	| Knstrn|
| Gm15650	| Dok7, Htra3|
| Gm28513	| Rgs16, Teddm2, Glul|
| Gm42705	| none|
| 8030451A03Rik	| Akna, Tmem268, Pappa|
| Gm10570	|Zbtb8os, Marcksl1, Fabp3, Zcchc17 |
| Gimap1os	| Rarres2, Aoc1, Igf2bp3|


## EdgeR

look at lncRNA genes across cell types and see how they cluster

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/heatmap%20top%20100%20lncrna%20all%20samples.png)

looking here there is reasonable clustering by cell type, suggesting a degree of cell type specificity of expression pattern. Gene names arent really legible as they are so small, but are uploaded as a separate list

### 1: Gm15527 (ENSMUSG00000086946)

Overlaps two CTCF bidning sites and multiple promoter regions

Overlapping protein coding genes: Chn2 (Decreased expression associated with gliomas and breast tumors,  increased expression is associated with lymphomas)

Human Homolog/syntenic location: CHN2, also overlaps with an antisense lncRNA

BLAST shows that this sequence is conserved in rodents and bats, but not in humans
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-20%20at%2016.00.47.png)

Behaviour of Gm15527 in different cell types

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |			
|	B	|	4	|	Gm15527	|	ENSMUSG00000086946	|	56.5843423	|	-3.132654655	|	0.861717019	|	-3.635363565	|	0.000277589	|	0.013534512	|	6	|	54222727	|	54256061	|	-1	|
|	LD	|	4	|	Gm15527	|	ENSMUSG00000086946	|	9.373327961	|	-4.338756071	|	1.381996683	|	-3.139483709	|	0.001692458	|	0.019286727	|	6	|	54222727	|	54256061	|	-1	|
|	LP	|	4	|	Gm15527	|	ENSMUSG00000086946	|	25.31616911	|	3.4950952	|	0.981346017	|	3.561531958	|	0.000368697	|	0.010260152	|		6	|	54222727	|	54256061	|	-1	|
|	S	|	4	|	Gm15527	|	ENSMUSG00000086946	|	46.33272281	|	-2.390907834	|	0.735387141	|	-3.251223335	|	0.001149096	|	0.021106651	|	6	|	54222727	|	54256061	|	-1	|

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/gm15527.png)

### 2: Gm15629 (ENSMUSG00000090054)

Overlapping protein coding genes: Lcp1 (ENSMUSG00000021998), expression is induced accompanying tumorigenesis in solid tissues. Rubcnl (ENSMUSG00000034959) also within 1mb (increased expression assocaited w/ cervical cancer, prognostic marker in endometrial cancer, strongly associated with gastric cancer https://www.proteinatlas.org/ENSG00000102445-RUBCNL/pathology)

Behaviour in Luminal differentiated cells - appears to be a nice correlation for all three - same patterns but much hihger count in lncRNA. According to ENCODE lncRNAs are generally expressed at a very low level relative to protein coding genes so this is interesting.

|		|		|	 	|
|	 ------------- |		 ------------- |		 ------------- |	
| ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2018.45.36.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2018.45.13.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2018.45.23.png) |
Human Homolog/syntenic location: Lcp1, associated w tumours in humans

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Mouse_Gm15629.png)

LncRNA has only been identified in two other Mus species

BLAST:
- first exon sequence appears to be unique to mice
- second exon in rodents generally
- third exon unique to mice
- BLAST of whole thing brings up Lcp1 as protein coding regions are conserved
 
Behaviour of  in different cell types:


|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		
|	B	|	4	|	Gm15629	|	ENSMUSG00000090054	|	256.9518598	|	-3.383026604	|	0.667130494	|	-5.071011798	|	3.96E-07	|	7.80E-05	|	14	|	75434513	|	75436711	|	-1	|
|	A	|	4	|	Gm15629	|	ENSMUSG00000090054	|	7.434596654	|	3.08077301	|	1.105515804	|	2.786729052	|	0.005324298	|	0.061857985		|	14	|	75434513	|	75436711	|	-1	|
|	LD	|	4	|	Gm15629	|	ENSMUSG00000090054	|	1079.308889	|	-2.472095177	|	0.37445647	|	-6.601822577	|	4.06E-11	|	7.81E-09		|	14	|	75434513	|	75436711	|	-1	|
|	S	|	4	|	Gm15629	|	ENSMUSG00000090054	|	573.9708165	|	-4.312986186	|	0.756863838	|	-5.698496837	|	1.21E-08	|	5.33E-06	|	14	|	75434513	|	75436711	|	-1	|

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/gm15629.png)

### 3: Gm42697 (ENSMUSG00000106120)

Overlapping protein coding genes: egf

Also overawlpping Egf: predicted gene 42650 (ENSMUSG00000105596) and 	RIKEN cDNA 6330410L21 gene (ENSMUSG00000105960)--> have a look at them too

also within 1 mb of this lncRNA (and signfificantly differnetially) expressed is Pla2g12a (ENSMUSG00000027999)

Human Homolog/syntenic location: lncRNA has orthologs in other mouse species, while Egf is conseved, the expression of an antisense lncRNA doesn't appear to be.

plotting expression in basal cells as this is the point that has the greatest padj and log2FoldChange
| | |
|	 ------------- |		 ------------- |
| ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2018.49.45.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2018.49.34.png) |
Behaviour in different cell types:

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |	
|	E	|	4	|	Gm42697	|	ENSMUSG00000106120	|	50.47332703	|	4.609328028	|	1.23496395	|	3.732358364	|	0.000189695	|	0.019549485	|	3	|	129522494	|	129533837	|	1	|
|	B	|	4	|	Gm42697	|	ENSMUSG00000106120	|	116.9573459	|	-17.21151338	|	1.605236214	|	-10.72210634	|	8.02E-27	|	1.42E-23		|	3	|	129522494	|	129533837	|	1	|
|	LD	|	4	|	Gm42697	|	ENSMUSG00000106120	|	88.48331024	|	-2.943625455	|	0.818927132	|	-3.594490084	|	0.000325028	|	0.005770489		|	3	|	129522494	|	129533837	|	1	|
|	LP	|	4	|	Gm42697	|	ENSMUSG00000106120	|	80.49656195	|	-2.880679672	|	0.870029301	|	-3.311014546	|	0.000929584	|	0.02237282		|	3	|	129522494	|	129533837	|	1	|

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/gm42697.png)


### 4: Gm13031 (ENSMUSG00000087698)

Overlapping protein coding genes: antisense to Padi2 (ENSMUSG00000028927)

Has four differentially expressed protein coding genes within 1mb
- Padi2 (ENSMUSG00000028927) developmentawlly important - severl cardiovascular morphology phenotypes, implicated in neurodegenerative stuff.
- Epha2 (ENSMUSG00000006445) abnormal development phenotype
- Fblim1 (ENSMUSG00000006219) abnormal development phenotypes
- Efhd2 (ENSMUSG00000040659)

| | | |
|	 ------------- |		 ------------- |		 ------------- |
|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2018.59.47.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2018.59.36.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2018.59.26.png) |
|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2018.59.17.png) | ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2018.59.07.png) | |

Looking at the expression here its interesting that the lncRNA starts very high and then crashes down whereas the others have little increases

Human Homolog/syntenic location: lncRNA conserved in mice. Protein coding gene conserved, but lncRNA not present
 
Behaviour in different cell types:

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |	
|	LD	|	3	|	Gm13031	|	ENSMUSG00000087698	|	94.31613416	|	5.131843056	|	0.649232894	|	7.904471729	|	2.69E-15	|	1.04E-12		|	4	|	140674977	|	140685213	|	-1	|
|	LP	|	3	|	Gm13031	|	ENSMUSG00000087698	|	55.68521887	|	3.456775209	|	0.76317455	|	4.529468665	|	5.91E-06	|	0.000421258	|	4	|	140674977	|	140685213	|	-1	|
|	S	|	3	|	Gm13031	|	ENSMUSG00000087698	|	132.1266243	|	-2.973004485	|	0.498134522	|	-5.96827635	|	2.40E-09	|	1.59E-06		|	4	|	140674977	|	140685213	|	-1	|

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/gm1301.png)

### 5: Gm14091 (ENSMUSG00000086421)

Overlapping protein coding genes: antisense to Knstrn ENSMUSG00000027331, which is mutated in 19% of cuteneous squamous cell carcinomas (Lee CS et al, Nat Genet 2014).

Human Homolog/syntenic location: Knstrn. but no antisense overlap lncRNA.

|||
|	 ------------- |		 ------------- |	
| ![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2019.07.59.png)|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2019.07.50.png)|
|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2019.12.44.png)|![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Screenshot%202021-02-21%20at%2019.12.36.png)|

looks almost like "opposite" expression patterns in stromal cells, but very similar patterns in luminal differentiated cells

cool

Behaviour in different cell types:

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |	
|	E	|	3	|	Gm14091	|	ENSMUSG00000086421	|	22.66956384	|	-2.050235133	|	0.732271129	|	-2.799830625	|	0.005112943	|	0.099562213		|	2	|	118660355	|	118664603	|	-1	|
|	LD	|	3	|	Gm14091	|	ENSMUSG00000086421	|	53.12878744	|	-2.398914535	|	0.508963665	|	-4.713331619	|	2.44E-06	|	0.000108165|	2	|	118660355	|	118664603	|	-1	|
|	S	|	3	|	Gm14091	|	ENSMUSG00000086421	|	24.62240463	|	-3.309827299	|	0.591046753	|	-5.599941594	|	2.14E-08	|	8.10E-06		|	2	|	118660355	|	118664603	|	-1	|

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/gm14091.png)


### 6: Gm15650 (ENSMUSG00000087703)

Overlapping protein coding genes: antisense to Sh3tc1 ENSMUSG00000036553

Human Homolog/syntenic location: protien coding gene conserved, lncRNA not

Behaviour in different cell types:

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |	
|	B	|	3	|	Gm15650	|	ENSMUSG00000087703	|	129.9594979	|	-2.437957064	|	0.592316328	|	-4.115971398	|	3.86E-05	|	0.002849873		|	5	|	35879899	|	35881829	|	1	|
|	LD	|	3	|	Gm15650	|	ENSMUSG00000087703	|	109.3269826	|	-2.338897108	|	0.812108942	|	-2.880028757	|	0.003976389	|	0.035710139|	5	|	35879899	|	35881829	|	1	|
|	S	|	3	|	Gm15650	|	ENSMUSG00000087703	|	67.395382	|	-2.932925134	|	0.847256959	|	-3.461671341	|	0.000536832	|	0.012565676		|	5	|	35879899	|	35881829	|	1	|

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/gm15650.png)


### 7: Gm28513 (ENSMUSG00000099568)

Overlapping protein coding genes: 
- antisense to Rgs16 ENSMUSG00000026475, which also has a sense lncRNA	predicted gene 28512 ENSMUSG00000100183
- antisense to Rnasel ENSMUSG00000066800

overlaps lots of promoter regions

Human Homolog/syntenic location: same protein, no lncRNA
 
Behaviourin different cell types:

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |	
|	B	|	3	|	Gm28513	|	ENSMUSG00000099568	|	297.8979686	|	-2.304077761	|	0.485387228	|	-4.746885845	|	2.07E-06	|	0.000229037	|	1	|	153619298	|	153625070	|	-1	|
|	LD	|	3	|	Gm28513	|	ENSMUSG00000099568	|	154.0247959	|	-3.789380639	|	1.144818015	|	-3.310028833	|	0.000932864	|	0.012303138	|	1	|	153619298	|	153625070	|	-1	|
|	LP	|	3	|	Gm28513	|	ENSMUSG00000099568	|	451.4280076	|	-5.255830861	|	0.605345175	|	-8.682370123	|	3.88E-18	|	6.90E-15		|	1	|	153619298	|	153625070	|	-1	|

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/gm28513.png)


### 8: Gm42705 (ENSMUSG00000105324)

Overlapping protein coding genes:no protein coding genes nearby, lots of other lncRNAs though. not conserved outside mice
http://rna.tbi.univie.ac.at//cgi-bin/RNAWebSuite/RNAfold.cgi?PAGE=3&ID=MXq_l_Hdki
had a look at folding

Behaviour in different cell types:

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |	
|	B	|	3	|	Gm42705	|	ENSMUSG00000105324	|	18.8978368	|	-2.872840932	|	1.018317909	|	-2.821163121	|	0.004784986	|	0.099865483		|	3	|	143449138	|	143532773	|	1	|
|	LD	|	3	|	Gm42705	|	ENSMUSG00000105324	|	6.209372551	|	-7.023785455	|	2.058735015	|	-3.411699614	|	0.000645592	|	0.009500153		|	3	|	143449138	|	143532773	|	1	|
|	LP	|	3	|	Gm42705	|	ENSMUSG00000105324	|	20.51994296	|	-2.442031993	|	0.873848629	|	-2.79457095	|	0.005196864	|	0.074513722		|	3	|	143449138	|	143532773	|	1	|

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/gm42705.png)


### 9: 8030451A03Rik (ENSMUSG00000073821)
http://rna.tbi.univie.ac.at//cgi-bin/RNAWebSuite/RNAfold.cgi?PAGE=3&ID=26D7JUEc6S

Overlapping protein coding genes:
- tenascin C ENSMUSG00000028364 (role in cancer)
- predicted gene 11216 ENSMUSG00000086755
- RIKEN cDNA 8030451A03 gene ENSMUSG00000073821

Human Homolog/syntenic location: human version of tnc also has a lncRNA, whose deletion is associated with cancer: DELEC1 ENSG00000173077 http://rna.tbi.univie.ac.at//cgi-bin/RNAWebSuite/RNAfold.cgi?PAGE=3&ID=4pHDjm2Ozq
BLAST 
Behaviour of  in different cell types:

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |	
|	B	|	3	|	8030451A03Rik	|	ENSMUSG00000073821	|	59.52842217	|	2.685963084	|	0.75757131	|	3.545492086	|	0.000391881	|	0.01782555		|	4	|	63898091	|	64068161	|	1	|
|	LD	|	3	|	8030451A03Rik	|	ENSMUSG00000073821	|	51.86269352	|	-3.016348472	|	1.065416227	|	-2.83114561	|	0.00463816	|	0.040061746	|	4	|	63898091	|	64068161	|	1	|
|	LP	|	3	|	8030451A03Rik	|	ENSMUSG00000073821	|	227.7021443	|	-3.064982509	|	0.502046385	|	-6.104978749	|	1.03E-09	|	3.66E-07		|	4	|	63898091	|	64068161	|	1	|

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/803...rik.png)

### 10: Gm10570 (ENSMUSG00000073752)

Overlapping protein coding genes:

Human Homolog/syntenic location: 
BLAST 
Behaviour of  in different cell types:

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |	
|	LD	|	3	|	Gm10570	|	ENSMUSG00000073752	|	7.831840165	|	-5.223209164	|	1.545952254	|	-3.378635498	|	0.000728465	|	0.010251816	|	4	|	130200349	|	130202467	|	-1	|
|	LP	|	3	|	Gm10570	|	ENSMUSG00000073752	|	17.76069584	|	-4.045936226	|	1.032762378	|	-3.917586767	|	8.94E-05	|	0.003250865	|	4	|	130200349	|	130202467	|	-1	|
|	S	|	3	|	Gm10570	|	ENSMUSG00000073752	|	3.619822916	|	3.097310447	|	1.162903431	|	2.663428762	|	0.007734879	|	0.081974977		|	4	|	130200349	|	130202467	|	-1	|


![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/gm10570.png)


### 11: Gimap1os (ENSMUSG00000044867)

Overlapping protein coding genes:

Human Homolog/syntenic location: 
BLAST 
Behaviour of  in different cell types:

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |	
|	LD	|	3	|	Gimap1os	|	ENSMUSG00000044867	|	40.66463525	|	-4.233494504	|	1.102813041	|	-3.838814327	|	0.00012363	|	0.002770269		|	6	|	48715280	|	48718378	|	-1	|
|	LP	|	3	|	Gimap1os	|	ENSMUSG00000044867	|	198.9553236	|	3.231920839	|	1.158399366	|	2.789988439	|	0.005270992	|	0.074513722		|	6	|	48715280	|	48718378	|	-1	|
|	S	|	3	|	Gimap1os	|	ENSMUSG00000044867	|	168.0123291	|	-3.90957232	|	0.71819819	|	-5.443584202	|	5.22E-08	|	1.38E-05		|	6	|	48715280	|	48718378	|	-1	|

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/gimap1os.png)
