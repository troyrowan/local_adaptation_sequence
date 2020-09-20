# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("rtracklayer")

library(tidyverse)
library(GALLO)

genes =
  import_gff_gtf(db_file = "data/Bos_taurus.ARS-UCD1.2.101.gtf.gz", file_type = "gtf")



gene_annotation = function(gwas_results, sig_filter = 1e-5, search_bp = 50000){
  genes =
    gwas_results %>%
      filter(p < sig_filter) %>%
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


read_gwas = function(gwas_path)
  ran_sex_gwas =
  read_tsv(gwas_path) %>%
  mutate(chrbp = paste(Chr, bp, sep = ":"),
         maf = case_when(Freq > 0.5 ~ 1-Freq,
                         TRUE ~ Freq))
