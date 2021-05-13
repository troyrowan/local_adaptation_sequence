library(furrr)
library(GenWin)
library(tidyverse)
#pacman::p_load(dplyr, furrr, tidyr, readr, GenWin)

genwin = function(selscandf, chr){
  tmpdf = filter(selscandf,
                 CHR == chr)
  windows =
    GenWin::splineAnalyze(Y = tmpdf$Mu,
                          map = tmpdf$BP) %>%
    .$windowData
  return(windows)
}


plan(multicore, workers = 1)

sim_raisd =
  seq(1,29) %>%
    future_map(
      ~read_delim(paste0("output/200907_SIM/sfs_selection/RAiSD_Report.200907_SIM.chr", .x,".gz"),
                 skip = 1,
                 delim = "\t",
                 col_names = c("BP", "stat", "finish", "VAR", "SFS", "LD", "Mu")) %>%
        mutate(CHR = .x,
               SNP = paste(CHR, BP, sep = ":"))
      ) %>%
    reduce(bind_rows)
head(sim_raisd)

ran_raisd =
  seq(1,29) %>%
    future_map(
      ~read_delim(paste0("output/200910_RAN/sfs_selection/RAiSD_Report.200910_RAN.chr", .x,".gz"),
                 skip = 1,
                 delim = "\t",
                 col_names = c("BP", "start", "finish", "VAR", "SFS", "LD", "Mu")) %>%
        mutate(CHR = .x,
               SNP = paste(CHR, BP, sep = ":"))) %>%
    reduce(bind_rows)
head(ran_raisd)

print("Running tester")
raisd_chr23 = read_delim("output/200907_SIM/sfs_selection/RAiSD_Report.200907_SIM.chr23.gz",
               skip = 1,
               delim = "\t",
               col_names = c("BP", "stat", "finish", "VAR", "SFS", "LD", "Mu")) %>%
  top_n(1000, wt = -BP) %>%
	mutate(CHR = 23)



genwin(raisd_chr23, 23)

print("calculating genome-wide windows for Red Angus")
# seq(1,29) %>%
#   future_map(
#     ~genwin(ran_raisd, .x) %>%
#       mutate(CHR = .x,
#              BP = (WindowStart + WindowStop)/2,
#              SNP = paste(CHR, BP, sep = ":")) %>%
#       select(SNP,CHR, BP, everything())
#   ) %>%
#   reduce(bind_rows) %>%
# 	write_tsv("output/200907_SIM/selscan_window/200907_SIM.RAiSD.windows.txt")
#
# seq(1,29) %>%
#   future_map(
#     ~genwin(ran_raisd, .x) %>%
#       mutate(CHR = .x,
#              BP = (WindowStart + WindowStop)/2,
#              SNP = paste(CHR, BP, sep = ":")) %>%
#       select(SNP,CHR, BP, everything())
#   ) %>%
#   reduce(bind_rows) %>%
# 	write_tsv("output/200910_RAN/selscan_window/200910_RAN.RAiSD.windows.txt")
