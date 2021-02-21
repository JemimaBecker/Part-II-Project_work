# Week: Identifying the best cell type and candidate genes

I ran DESeq2 on the lncRNA genes in each cell type individually to examine:

- How many genes were differentially expressed
- What the distribution of expression looked like

| Cell type  |Number of lncRNA genes with padj<0.1 and log2FoldChange > 2|
| ------------- | ------------- |
|Luminal Differentiated| 234 |
|Stromal|132|
|Adipocyte|66|
|Luminal Progenitor|62|
|Basal|40|
|Endothelial|23|

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
| Gm15629	| 3| - |- | + | + | +|- |	
| Gm42697	| 4|  +| + |-  | + | -| + |	
| Gm13031	| 3| + | + | + |- |- |- |	
| Gm14091	| 3| + |- | + |- | -|  +|	
| Gm15650	| 3| + |- |- | + |- |- |	
| Gm28513	| 3| + | + |- | -|- |- |	
| Gm42705	| 3| + | + |- | + | -|- |	
| 8030451A03Rik	| 3| + | + |- | + | |- |	
| Gm10570	| 3| + | + |- | -| -|- |	
| Gimap1os	| 3| + | + |- | -|- | -|	


Time to look at the behaviour of each of these lncRNAs in their respective cell types 

### 1: Gm15527 (ENSMUSG00000086946)

Overlaps two CTCF bidning sites and multiple promoter regions

Overlapping protein coding genes: Chn2 (Decreased expression associated with gliomas and breast tumors,  increased expression is associated with lymphomas)

Human Homolog/syntenic location: CHN2, also overlaps with an antisense lncRNA

BLAST shows that this sequence is conserved in rodents and bats, but not in humans
![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/images/Screenshot%202021-02-20%20at%2016.00.47.png)

Behaviour of Gm15527 in different cell types

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	description	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |	
|	B	|	4	|	Gm15527	|	ENSMUSG00000086946	|	56.5843423	|	-3.132654655	|	0.861717019	|	-3.635363565	|	0.000277589	|	0.013534512	|	predicted gene 15527 [Source:MGI Symbol;Acc:MGI:3782974]	|	6	|	54222727	|	54256061	|	-1	|
|	LD	|	4	|	Gm15527	|	ENSMUSG00000086946	|	9.373327961	|	-4.338756071	|	1.381996683	|	-3.139483709	|	0.001692458	|	0.019286727	|	predicted gene 15527 [Source:MGI Symbol;Acc:MGI:3782974]	|	6	|	54222727	|	54256061	|	-1	|
|	LP	|	4	|	Gm15527	|	ENSMUSG00000086946	|	25.31616911	|	3.4950952	|	0.981346017	|	3.561531958	|	0.000368697	|	0.010260152	|	predicted gene 15527 [Source:MGI Symbol;Acc:MGI:3782974]	|	6	|	54222727	|	54256061	|	-1	|
|	S	|	4	|	Gm15527	|	ENSMUSG00000086946	|	46.33272281	|	-2.390907834	|	0.735387141	|	-3.251223335	|	0.001149096	|	0.021106651	|	predicted gene 15527 [Source:MGI Symbol;Acc:MGI:3782974]	|	6	|	54222727	|	54256061	|	-1	|

### 2: Gm15629 (ENSMUSG00000090054)

Overlapping protein coding genes: Lcp1 (ENSMUSG00000021998), expression is induced accompanying tumorigenesis in solid tissues.
Human Homolog/syntenic location: Lcp1, associated w tumours in humans

![](https://github.com/AFS-Part-II-Projects/Jemima_Becker/blob/main/Week%205%20Images/Mouse_Gm15629.png)

LncRNA has only been identified in two other Mus species
BLAST ?? cant find in genbank for some reason
Behaviour of  in different cell types:

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	description	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |	
|	B	|	4	|	Gm15629	|	ENSMUSG00000090054	|	256.9518598	|	-3.383026604	|	0.667130494	|	-5.071011798	|	3.96E-07	|	7.80E-05	|	predicted gene 15629 [Source:MGI Symbol;Acc:MGI:3783073]	|	14	|	75434513	|	75436711	|	-1	|
|	A	|	4	|	Gm15629	|	ENSMUSG00000090054	|	7.434596654	|	3.08077301	|	1.105515804	|	2.786729052	|	0.005324298	|	0.061857985	|	predicted gene 15629 [Source:MGI Symbol;Acc:MGI:3783073]	|	14	|	75434513	|	75436711	|	-1	|
|	LD	|	4	|	Gm15629	|	ENSMUSG00000090054	|	1079.308889	|	-2.472095177	|	0.37445647	|	-6.601822577	|	4.06E-11	|	7.81E-09	|	predicted gene 15629 [Source:MGI Symbol;Acc:MGI:3783073]	|	14	|	75434513	|	75436711	|	-1	|
|	S	|	4	|	Gm15629	|	ENSMUSG00000090054	|	573.9708165	|	-4.312986186	|	0.756863838	|	-5.698496837	|	1.21E-08	|	5.33E-06	|	predicted gene 15629 [Source:MGI Symbol;Acc:MGI:3783073]	|	14	|	75434513	|	75436711	|	-1	|

### 3: Gm42697 (ENSMUSG00000106120)

Overlapping protein coding genes:

Human Homolog/syntenic location: 
BLAST 
Behaviour of  in different cell types:

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	description	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |	
|	E	|	4	|	Gm42697	|	ENSMUSG00000106120	|	50.47332703	|	4.609328028	|	1.23496395	|	3.732358364	|	0.000189695	|	0.019549485	|	predicted gene 42697 [Source:MGI Symbol;Acc:MGI:5662834]	|	3	|	129522494	|	129533837	|	1	|
|	B	|	4	|	Gm42697	|	ENSMUSG00000106120	|	116.9573459	|	-17.21151338	|	1.605236214	|	-10.72210634	|	8.02E-27	|	1.42E-23	|	predicted gene 42697 [Source:MGI Symbol;Acc:MGI:5662834]	|	3	|	129522494	|	129533837	|	1	|
|	LD	|	4	|	Gm42697	|	ENSMUSG00000106120	|	88.48331024	|	-2.943625455	|	0.818927132	|	-3.594490084	|	0.000325028	|	0.005770489	|	predicted gene 42697 [Source:MGI Symbol;Acc:MGI:5662834]	|	3	|	129522494	|	129533837	|	1	|
|	LP	|	4	|	Gm42697	|	ENSMUSG00000106120	|	80.49656195	|	-2.880679672	|	0.870029301	|	-3.311014546	|	0.000929584	|	0.02237282	|	predicted gene 42697 [Source:MGI Symbol;Acc:MGI:5662834]	|	3	|	129522494	|	129533837	|	1	|


### 4: Gm13031 (ENSMUSG00000087698)

Overlapping protein coding genes:

Human Homolog/syntenic location: 
BLAST 
Behaviour of  in different cell types:

|	Cell type	|	Number of cell types	|	external_gene_name	|	ensembl_gene_id	|	baseMean	|	log2FoldChange	|	lfcSE	|	stat	|	pvalue	|	padj	|	description	|	chromosome_name	|	start_position	|	end_position	|	strand	|
|	 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |		 ------------- |	
|	LD	|	3	|	Gm13031	|	ENSMUSG00000087698	|	94.31613416	|	5.131843056	|	0.649232894	|	7.904471729	|	2.69E-15	|	1.04E-12	|	predicted gene 13031 [Source:MGI Symbol;Acc:MGI:3702680]	|	4	|	140674977	|	140685213	|	-1	|
|	LP	|	3	|	Gm13031	|	ENSMUSG00000087698	|	55.68521887	|	3.456775209	|	0.76317455	|	4.529468665	|	5.91E-06	|	0.000421258	|	predicted gene 13031 [Source:MGI Symbol;Acc:MGI:3702680]	|	4	|	140674977	|	140685213	|	-1	|
|	S	|	3	|	Gm13031	|	ENSMUSG00000087698	|	132.1266243	|	-2.973004485	|	0.498134522	|	-5.96827635	|	2.40E-09	|	1.59E-06	|	predicted gene 13031 [Source:MGI Symbol;Acc:MGI:3702680]	|	4	|	140674977	|	140685213	|	-1	|

### 5: Gm14091 ()

Overlapping protein coding genes:

Human Homolog/syntenic location: 
BLAST 
Behaviour of  in different cell types:

| Cell type  | log2FoldChange| Adjusted P-value |
| ------------- | ------------- |------------- |
|Luminal Differentiated||
|Luminal Progenitor|||
|Basal|||
|Stromal|||

### 6: Gm15650 ()

Overlapping protein coding genes:

Human Homolog/syntenic location: 
BLAST 
Behaviour of  in different cell types:

| Cell type  | log2FoldChange| Adjusted P-value |
| ------------- | ------------- |------------- |
|Luminal Differentiated||
|Luminal Progenitor|||
|Basal|||
|Stromal|||

### 7: Gm28513 ()

Overlapping protein coding genes:

Human Homolog/syntenic location: 
BLAST 
Behaviour of  in different cell types:

| Cell type  | log2FoldChange| Adjusted P-value |
| ------------- | ------------- |------------- |
|Luminal Differentiated||
|Luminal Progenitor|||
|Basal|||
|Stromal|||

### 8: Gm42705 ()

Overlapping protein coding genes:

Human Homolog/syntenic location: 
BLAST 
Behaviour of  in different cell types:

| Cell type  | log2FoldChange| Adjusted P-value |
| ------------- | ------------- |------------- |
|Luminal Differentiated||
|Luminal Progenitor|||
|Basal|||
|Stromal|||

### 9: 8030451A03Rik ()

Overlapping protein coding genes:

Human Homolog/syntenic location: 
BLAST 
Behaviour of  in different cell types:

| Cell type  | log2FoldChange| Adjusted P-value |
| ------------- | ------------- |------------- |
|Luminal Differentiated||
|Luminal Progenitor|||
|Basal|||
|Stromal|||

### 10: Gm10570 ()

Overlapping protein coding genes:

Human Homolog/syntenic location: 
BLAST 
Behaviour of  in different cell types:

| Cell type  | log2FoldChange| Adjusted P-value |
| ------------- | ------------- |------------- |
|Luminal Differentiated||
|Luminal Progenitor|||
|Basal|||
|Stromal|||

### 11: Gimap1os ()

Overlapping protein coding genes:

Human Homolog/syntenic location: 
BLAST 
Behaviour of  in different cell types:

| Cell type  | log2FoldChange| Adjusted P-value |
| ------------- | ------------- |------------- |
|Luminal Differentiated||
|Luminal Progenitor|||
|Basal|||
|Stromal|||

