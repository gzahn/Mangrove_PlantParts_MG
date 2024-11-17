# SETUP ####

## packages ####
library(tidyverse)
library(janitor)
library(DESeq2) # https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

## data ####
mat <- readRDS("./data/count_table.RDS")
sample_meta <- readRDS("./data/sample_meta.RDS")
uniprot_meta <-readRDS("./data/gene_meta.RDS")


# BUILD DESEQ ####

# build DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = sample_meta,
                              design = ~ plant_part)
dds

# add gene metadata to deseq object
mcols(dds) <- DataFrame(mcols(dds),uniprot_meta)
dds

# Differential Expression analysis
dds <- DESeq(dds)
res <- results(dds,alpha = 0.05)
resultsNames(dds)
res

# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="plant_part_Leaf_vs_Sediment", type="apeglm")

# quick viz
plotMA(resLFC)

# pull significant genes from leaf vs sediment
sig_genes_leaf_vs_sediment <- 
  resLFC %>%
  as.data.frame() %>% 
  dplyr::filter(padj < 0.01 & abs(log2FoldChange) > 2)
sig_genes_leaf_vs_sediment$uniprot_id <- row.names(sig_genes_leaf_vs_sediment)

# combine with gene info
sig_genes_leaf_vs_sediment <- 
  uniprot_meta %>% 
  dplyr::filter(uniprot_id %in% row.names(sig_genes_leaf_vs_sediment)) %>% 
  full_join(sig_genes_leaf_vs_sediment)

sig_genes_leaf_vs_sediment %>%
  dplyr::filter(log2FoldChange > 0 & !is.na(pfam_id)) %>% # just look at positive values for now (enriched in leaf)
  ggplot(aes(x=ko_description,y=log2FoldChange)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust = .5))

# need to subset to something of interest, like CAZy genes or a list of ko genes in C cycle
