#SETUP ####
library(tidyverse)
library(janitor)
library(DESeq2) # https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


# functions
read_stats <- function(x){
  sampleid <- x %>% basename() %>% str_remove("_scafstats.txt|_refstats.txt")
  dat <- read_delim(x) %>% clean_names() %>% mutate(sample_id=sampleid)
  names(dat)[1] <- "full_name"
  dat <- 
    dat %>% 
    mutate(uniprot_id = full_name %>% str_split("\\|") %>% map_chr(1),
           predicted_gene = full_name %>% str_split("\\|") %>% map_chr(2)) %>% 
    dplyr::select(sample_id, predicted_gene, uniprot_id, unambiguous_reads)
  return(dat)
}

read_annot <- function(x){
  dat <- read_delim(x,delim = "\t",skip = 1,col_names = FALSE)
  sampleid <- x %>% basename() %>% str_remove("_annotation_table.tsv")
  header <- readLines(x)[1] %>% 
    str_squish() %>% 
    str_split("\\s") %>% 
    unlist()
  colnames(dat) <- header
  dat$sample_id=sampleid
  return(dat)
}

# file paths
path <- "./data/annotations"
# scaffold stats
scafstats_fp <- list.files(path,full.names = TRUE, pattern = "scafstats.txt")
# reference stats (overall for that file, might be handy)
refstats_fp <- list.files(path,full.names = TRUE, pattern = "refstats.txt")
# annotation tables
annot_fp <- list.files(path,full.names = TRUE, pattern = "annotation_table.tsv")
# sample names
sample_names <- basename(annot_fp) %>% str_remove("_annotation_table.tsv")
# subset to remove scafstats files that have no corresponding annotation file
scafstats_fp <- grep(x = scafstats_fp,pattern = paste(sample_names,collapse = "|"),value = TRUE)


# read in delim files
scafstats <- map(scafstats_fp,read_stats)
refstats <- map(refstats_fp,read_delim) %>% purrr::reduce(full_join) %>% clean_names()
names(scafstats) <- sample_names

annotations <- map(annot_fp,read_annot)

# stick them all together
dat <- 
  full_join(
    purrr::reduce(scafstats,full_join),
    purrr::reduce(annotations,full_join)
  )

 head(dat)

# build gene_metadata
uniprot_meta <- 
  dat %>% 
  dplyr::select(uniprot_id,ko_id,ko_description,pfam_id,pfam_description,go_term) %>% 
  unique.data.frame()


# build sample metadata
sample_meta <- 
data.frame(
  sample_id = dat$sample_id %>% unique
) %>% 
  mutate(location = sample_id %>% 
           str_split("_") %>% 
           map_chr(1),
         location = case_when(location == "CJ" ~ "Chek Jawa",
                              location == "Sem" ~ "Semaku"),
         plant_part = sample_id %>% 
           str_split("_") %>% map_chr(2),
         plant_part = case_when(plant_part == "L" ~ "Leaf",
                                plant_part == "Pn" ~ "Pneumatophore",
                                plant_part == "Sed" ~ "Sediment"),
         replicate = sample_id %>% 
           str_split("_") %>% map_chr(3))
# enforce factor order with sediment as intercept
sample_meta$plant_part <- factor(sample_meta$plant_part,levels = c("Sediment","Pneumatophore","Leaf"))


# build count table (full)
count_tab <- 
  dat %>% 
  dplyr::select(sample_id, uniprot_id, unambiguous_reads) %>% 
  distinct(sample_id,uniprot_id, .keep_all = TRUE) %>% 
  pivot_wider(
              names_from = sample_id,
              values_from = unambiguous_reads)
# remove uniprot_ids that had no unambiguous mappings?
count_tab <- 
  count_tab[count_tab %>% 
            dplyr::select(-uniprot_id) %>% 
            rowSums(na.rm = TRUE) > 0,]


# convert to matrix
mat <- 
  count_tab %>% 
  dplyr::select(-uniprot_id) %>% 
  as.matrix
row.names(mat) <- count_tab$uniprot_id

# Convert NA to 0
mat[is.na(mat)] <- 0

# make sure gene metadata matches rows of count table (mat)
uniprot_meta_distinct <- distinct(uniprot_meta,uniprot_id,.keep_all = TRUE)
rownames(uniprot_meta_distinct) <- uniprot_meta_distinct$uniprot_id
uniprot_meta_distinct <- uniprot_meta_distinct[rownames(mat),]

# check order of input data
if(
  identical(
    colnames(mat),
    sample_meta$sample_id
  ) 
  &
  identical(
    rownames(mat),
    uniprot_meta_distinct$uniprot_id
  )
  
){
  # save count table, sample metadata, and gene metadata objects
  saveRDS(mat,"./data/count_table.RDS")
  saveRDS(sample_meta,"./data/sample_meta.RDS")
  saveRDS(uniprot_meta_distinct,"./data/gene_meta.RDS")
} else {
  cat("check that sample and gene metadata match the gene count matrix.")
}






# # build DESeqDataSet
# dds <- DESeqDataSetFromMatrix(countData = mat,
#                               colData = sample_meta,
#                               design = ~ plant_part)
# dds
# # add gene feature data
# uniprot_meta_distinct <- distinct(uniprot_meta,uniprot_id,.keep_all = TRUE)
# identical(uniprot_meta_distinct$uniprot_id,rownames(dds))
# mcols(dds) <- DataFrame(mcols(dds),uniprot_meta_distinct)
# dds
# 
# # Differential Expression analysis
# dds <- DESeq(dds)
# res <- results(dds,alpha = 0.05)
# resultsNames(dds)
# res
# 
# # Log fold change shrinkage for visualization and ranking
# resLFC <- lfcShrink(dds, coef="plant_part_Sediment_vs_Leaf", type="apeglm")
# 
# # quick viz
# plotMA(res)
