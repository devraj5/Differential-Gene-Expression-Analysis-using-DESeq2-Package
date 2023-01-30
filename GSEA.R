library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(DOSE)
library(dplyr)
library(tidyverse)

# SET THE DESIRED ORGANISM HERE
organism = "org.Mm.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)


# reading in data from deseq2
df = read.csv(file.choose(), header=TRUE)



# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$X

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
view(gse)

dotplot(gse, showCategory=10, split=".sign") + 
  ggtitle("Gene Set enrichment analysis") +
  theme(axis.text.y = element_text(size = rel(1))) +
  facet_grid(.~.sign)

gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)