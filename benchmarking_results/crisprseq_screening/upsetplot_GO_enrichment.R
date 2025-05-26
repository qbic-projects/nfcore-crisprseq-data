library(BiocManager)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(UpSetR)


table_bagel2 <- read.delim("./crisprseq_screening_results/crisprseq_test_full_wo_crisprcleanr_hitselection.1/bagel2/precision_recall/Brunello_RepA_Dropout_A375,Brunello_RepB_Dropout_A375_vs_Brunello_pDNA.tsv",
           header = TRUE, sep = "\t")
table_mageck_crisprseq <- read.delim("./crisprseq_screening_results/crisprseq_test_full_wo_crisprcleanr_hitselection.1/mageck/mle/Brunello_RepA_Dropout_A375,Brunello_RepB_Dropout_A375_vs_Brunello_pDNA/Brunello_RepA_Dropout_A375,Brunello_RepB_Dropout_A375_vs_Brunello_pDNA.gene_summary.txt",
                                     header = TRUE, sep = "\t")
table_mageck_vispr <- read.delim("mageck_vispr_results/results/results_vispr/test/mle.gene_summary.txt")


filtered_mageck <- table_mageck_crisprseq[table_mageck_crisprseq$Brunello_RepA_Dropout_A375_Brunello_RepB_Dropout_A375_vs_Brunello_pDNA.fdr
                                          < 0.1, ]
genes_mageck <- filtered_mageck$Gene


filtered_bagel2 <- table_bagel2[table_bagel2$FDR < 0.1, ]
genes_bagel2 <- filtered_bagel2$Gene 

exclusive_bagel2_genes <- setdiff(genes_bagel2, genes_mageck)


####Gene enrichment
# Convert Gene Symbols to ENTREZ IDs
gene_conversion <- bitr(exclusive_bagel2_genes, 
                        fromType = "SYMBOL",  # Input: Gene Symbols
                        toType = "ENTREZID",  # Output: ENTREZ IDs
                        OrgDb = org.Hs.eg.db)


ego <- enrichGO(gene = gene_conversion$ENTREZID, 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENTREZID",
                ont = "BP", 
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)
# Dot plot visualization
dotplot(ego, showCategory = 10)

# Mageck genes

exclusive_mageck_genes <- setdiff(genes_mageck, genes_bagel2)
gene_conversion_m <- bitr(exclusive_mageck_genes, 
                        fromType = "SYMBOL",  # Input: Gene Symbols
                        toType = "ENTREZID",  # Output: ENTREZ IDs
                        OrgDb = org.Hs.eg.db)
ego_m <- enrichGO(gene = gene_conversion_m$ENTREZID, 
                OrgDb = org.Hs.eg.db, 
                keyType = "ENTREZID",
                ont = "BP", 
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)
dotplot(ego_m, showCategory = 10)



######Upset plot
genes_vispr <- table_mageck_vispr[table_mageck_vispr$A375_SKIN.fdr < 0.1, ]$Gene

gene_lists <- list(
       "MAGeCK-VISPR MAGeCK MLE" = genes_vispr,
       "nf-core/crisprseq MAGeCK MLE"  = genes_mageck,
       "nf-core/crisprseq BAGEL2" = genes_bagel2
   )

upset(fromList(gene_lists), order.by = "freq",nsets = 4)

