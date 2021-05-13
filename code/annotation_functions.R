# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# BiocManager::install("rtracklayer")

library(tidyverse)
library(GALLO)
setwd(here::here())
genes =
  import_gff_gtf(db_file = "data/Bos_taurus.ARS-UCD1.2.103.gtf.gz", file_type = "gtf")

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


genwin = function(selscandf, chr){
  tmpdf = filter(selscandf,
                 CHR == chr)
  windows =
    GenWin::splineAnalyze(Y = tmpdf$norm_nsl,
                          map = tmpdf$BP) %>%
    .$windowData
  return(windows)
}


windowfinder = function(chromosome, basepair, windows){
  window_id =
    windows %>%
    filter(CHR == chromosome &
             WindowStart <= basepair &
             WindowStop >= basepair) %>%
    .$SNP
  return(window_id)
}


# FAANG =
#   read_tsv("data/FAANG_peaks.bed.gz",
#            col_names = c("CHR", "START", "STOP", "TMARK", "SCORE", "A", "B", "C", "D")) %>%
#   mutate(CHR = as.numeric(str_replace(CHR, "chr",""))) %>%
#   mutate(tmark = str_split_fixed(TMARK, "_", n = 5)[,1],
#          MARK = str_replace(tmark, "Macs2/", ""),
#          TISSUE = str_split_fixed(TMARK, "_", n = 5)[,2],
#          REP = str_split_fixed(TMARK, "_", n = 5)[,3]) %>%
#   select(CHR, START, STOP, MARK, TISSUE, REP)

# ATAC_peaks =
#   read_tsv("data/ATAC_peaks_ARS.bed",
#            col_names = c("CHR", "start", "stop", "name", "score", "strand", paste0("X", seq(1,8)))) %>%
#   mutate(CHR = as.numeric(str_replace(CHR, "chr", "")),
#          width = stop - start,
#          av_score = (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8)/8)
#
# FAANG_ATAC_peaks =
#   read_tsv("data/FAANG_ATAC_peaks_ARS.bed",
#            col_names = c("CHR", "start", "stop", "name", "score", "strand", paste0("X", seq(1,8)))) %>%
#   mutate(CHR = as.numeric(str_replace(CHR, "chr", "")),
#          width = stop - start,
#          av_score = (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8)/8)
#
# FAETH_Scores =
#   vroom("data/FAETH.txt.gz",
#         delim = "\t") %>%
#   filter(!is.na(FAETH))
