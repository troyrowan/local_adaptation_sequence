library(dplyr)
library(furrr)
library(tidyr)
library(readr)
library(GenWin)
pacman::p_load(dplyr, furrr, tidyr, readr, GenWin)


genwin = function(selscandf, chr){
  tmpdf = filter(selscandf,
                 CHR == chr)
  windows =
    GenWin::splineAnalyze(Y = tmpdf$norm_nsl,
                          map = tmpdf$BP) %>%
    .$windowData
  return(windows)
}

n_cores = availableCores() - 1

plan(multicore, workers = 30)

pb_sim_nsl =
  read_tsv("200907_SIM.purebred.nsl.our.100bins.norm.gz",
           col_names = c("SNP", "BP", "freq", "sl1", "sl0", "nsl", "norm_nsl", "sig")) %>%
  mutate(CHR = as.numeric(str_split_fixed(SNP, ":", n = 4)[,1]),
         SNP = paste(CHR, BP, sep = ":"))


seq(1,29) %>%
  future_map(
    ~genwin(pb_sim_nsl, .x) %>%
      mutate(CHR = .x,
             BP = (WindowStart + WindowStop)/2,
             SNP = paste(CHR, BP, sep = ":")) %>%
      select(SNP,CHR, BP, everything())
  ) %>%
  reduce(bind_rows) %>%
  write_tsv("200907_SIM.purebred.nsl.windows.txt")



all_sim_nsl =
  read_tsv("200907_SIM.somesim.nsl.our.100bins.norm.gz",
           col_names = c("SNP", "BP", "freq", "sl1", "sl0", "nsl", "norm_nsl", "sig")) %>%
  mutate(CHR = as.numeric(str_split_fixed(SNP, ":", n = 4)[,1]),
         SNP = paste(CHR, BP, sep = ":"))


seq(1,29) %>%
  future_map(
    ~genwin(all_sim_nsl, .x) %>%
      mutate(CHR = .x,
             BP = (WindowStart + WindowStop)/2,
             SNP = paste(CHR, BP, sep = ":")) %>%
      select(SNP,CHR, BP, everything())
  ) %>%
  reduce(bind_rows) %>%
  write_tsv("200907_SIM.all.nsl.windows.txt")



ran_nsl =
  read_tsv("200910_RAN.all.nsl.our.100bins.norm.gz",
           col_names = c("SNP", "BP", "freq", "sl1", "sl0", "nsl", "norm_nsl", "sig")) %>%
  mutate(CHR = as.numeric(str_split_fixed(SNP, ":", n = 4)[,1]),
         SNP = paste(CHR, BP, sep = ":"))

seq(1,29) %>%
  future_map(
    ~genwin(ran_nsl, .x) %>%
      mutate(CHR = .x,
             BP = (WindowStart + WindowStop)/2,
             SNP = paste(CHR, BP, sep = ":")) %>%
      select(SNP,CHR, BP, everything())
  ) %>%
  reduce(bind_rows) %>%
  write_tsv("200910_RAN.all.nsl.windows.txt")
