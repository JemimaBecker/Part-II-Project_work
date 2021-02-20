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
