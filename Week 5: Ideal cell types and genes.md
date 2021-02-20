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
| Gm15629	| 4| + |- | + | + | -|- |	
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

| Cell type  | log2FoldChange| Adjusted P-value |
| ------------- | ------------- |------------- |
|Luminal Differentiated|-4.338756|0.01928673|
|Luminal Progenitor|3.495095|0.01026015|
|Basal|-3.132655|0.01353451|
|Stromal|-2.390908|0.02110665|

### 2: Gm15629 ()

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

### 3: Gm42697 ()

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

### 4: Gm13031 ()

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

