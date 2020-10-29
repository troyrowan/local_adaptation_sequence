# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("rtracklayer")

library(tidyverse)
library(GALLO)
setwd(here::here())
genes =
  import_gff_gtf(db_file = "data/Bos_taurus.ARS-UCD1.2.101.gtf.gz", file_type = "gtf")

qtl =
  import_gff_gtf(db_file = "data/Bos_taurus.ARS-UCD1.2.QTL.gff.gz", file_type = "gff")

gene_annotation = function(gwas_results, search_bp = 50000){
  genes =
    gwas_results %>%
      dplyr::select(CHR, BP) %>%
      find_genes_qtls_around_markers(db_file=genes,
                                     marker_file=., method = "gene",
                                     marker = "snp", interval = search_bp, nThreads = NULL) %>%
      group_by(CHR, BP, gene_id) %>%
      mutate(upstream = start_pos - BP,
             downstream = BP - end_pos,
             distance = case_when(BP > start_pos & BP < end_pos ~ 0,
                                  TRUE ~ min (abs(upstream), abs(downstream)))) %>%
      ungroup()
  return(genes)
}

qtl_annotation = function(gwas_results, search_bp = 50000){
  qtl =
    gwas_results %>%
    dplyr::select(CHR, BP) %>%
    find_genes_qtls_around_markers(db_file=qtl,
                                   marker_file=., method = "qtl",
                                   marker = "snp", interval = search_bp, nThreads = NULL)
  return(qtl)
}

read_bolt_gwas = function(gwas_path){
  gwas_results =
    read_tsv(gwas_path) %>%
    mutate(chrbp = paste(CHR, BP, sep = ":"),
           maf = case_when(A1FREQ > 0.5 ~ 1-A1FREQ,
                           TRUE ~ A1FREQ),
           p = P_BOLT_LMM_INF)
  return(gwas_results)
}
