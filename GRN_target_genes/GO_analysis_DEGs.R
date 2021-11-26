library(DOSE)
library(clusterProfiler)
library(org.Mm.eg.db)

# Load list of DEGs - change files base on plot 
# For positive markers of TC electroporated, load - TC_pos_markers_Elect.csv
# For DEGs in AP, load - AP_DEGs_Elect.csv
# For positive markers of TC microinjected, load - TC_pos_markers_aMI.csv
# Load list of genes decreased in expression in TCs vs APs, cutoff Padj < 0.2
symbols <- as.character(read.csv2("AP_DEGs_Elect.csv")[,1])#AP_DEGs_Elect.csv/DEGs_in_APs_decreased_in_EPZ.csv
#data("geneList", package = "DOSE")
gene <- mapIds(org.Mm.eg.db, symbols, 'ENTREZID', 'SYMBOL')

library(enrichplot)
barplot(ggo, showCategory=20)
edo <- enrichDGN(gene)
barplot(edo,showCategory = 20)
#edo2 <- gseNCG(gene, nPerm=10000)
#dotplot(ego, showCategory= 20,orderBy="GeneRatio",title = "downregulated genes in TCs")


ego <- enrichGO(gene          = gene,
                OrgDb         = org.Mm.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
dotplot(ego, showCategory= 20,orderBy="GeneRatio",title = "downregulated genes APs_EPZ - CC")
barplot(ego,showCategory = 20,title = "downregulated genes APs_EPZ - CC")

ego1 <- enrichGO(gene          = gene,
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego1)
dotplot(ego1, showCategory= 20,orderBy="GeneRatio",title = "DEGs APs_EPZ - BP")
barplot(ego1,showCategory = 20,title = "downregulated genes APs_EPZ - BP")


ego2 <- enrichGO(gene          = gene,
                 OrgDb         = org.Mm.eg.db,
                 ont           = "ALL",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.05,
                 readable      = TRUE)
head(ego2)
dotplot(ego2, showCategory= 20,orderBy="GeneRatio",title = "DEGs genes APs -EPZ vs DMSO")
barplot(ego2,showCategory = 20,title = "downregulated genes APs_EPZ - BP")









