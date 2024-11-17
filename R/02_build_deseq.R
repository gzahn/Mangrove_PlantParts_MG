# SETUP ####

## packages ####
library(tidyverse)
library(janitor)
library(DESeq2) # https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

## data ####
mat <- readRDS("./data/count_table.RDS")
sample_meta <- readRDS("./data/sample_meta.RDS")
gene_meta <-readRDS("./data/gene_meta.RDS")

# BUILD DESEQ ####

# build DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = mat,
                              colData = sample_meta,
                              design = ~ plant_part)
dds
# add gene feature data
uniprot_meta_distinct <- distinct(uniprot_meta,uniprot_id,.keep_all = TRUE)
identical(uniprot_meta_distinct$uniprot_id,rownames(dds))
mcols(dds) <- DataFrame(mcols(dds),uniprot_meta_distinct)
dds

# Differential Expression analysis
dds <- DESeq(dds)
res <- results(dds,alpha = 0.05)
resultsNames(dds)
res

# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="plant_part_Sediment_vs_Leaf", type="apeglm")

# quick viz
plotMA(res)
