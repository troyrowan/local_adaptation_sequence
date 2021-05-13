library(furrr)
library(GenWin)
library(tidyverse)
#pacman::p_load(dplyr, furrr, tidyr, readr, GenWin)


genwin = function(selscandf, chr){
  tmpdf = filter(selscandf,
                 CHR == chr)
  windows =
    GenWin::splineAnalyze(Y = abs(tmpdf$norm_xpehh),
                          map = tmpdf$BP) %>%
    .$windowData
  return(windows)
}

n_cores = availableCores() - 1

plan(multicore, workers = 30)


sim_xpnsl =
	read_tsv("../200907_SIM/selscan/xpnsl/200907_SIM.xpnsl.out.norm.gz",
					 col_names = c("SNP", "BP", "gpos", "p1", "sl1", "p2", "sl2", "xpnsl", "norm_xpnsl", "sig"),
					 col_type = cols(SNP = col_character()),
					 skip = 1) %>%
	  mutate(CHR = as.numeric(str_split_fixed(SNP, ":", n = 4)[,1]),
	         SNP = paste(CHR, BP, sep = ":"))

print(head(sim_xpnsl))

sim_xpnsl_windows =
	seq(1,29) %>%
	  future_map(
	    ~genwin(sim_xpnsl, .x) %>%
	      mutate(CHR = .x,
	             BP = (WindowStart + WindowStop)/2,
	             SNP = paste(CHR, BP, sep = ":")) %>%
	      select(SNP,CHR, BP, everything())
	  ) %>%
	  reduce(bind_rows)

sim_xpnsl %>%
  mutate(window = map2_chr(.x = CHR, .y = BP,
                           ~ windowfinder(.x, .y, sim_xpnsl_windows)[1])) %>%
	write_tsv("200907_SIM.xpnsl.out.norm.windows.txt")







	ran_xpnsl =
	read_tsv("../200910_RAN/selscan/xpnsl/200907_RAN.xpnsl.out.norm.gz",
					 col_names = c("SNP", "BP", "gpos", "p1", "sl1", "p2", "sl2", "xpnsl", "norm_xpnsl", "sig"),
					 col_type = cols(SNP = col_character()),
					 skip = 1) %>% 
	  mutate(CHR = as.numeric(str_split_fixed(SNP, ":", n = 4)[,1]),
	         SNP = paste(CHR, BP, sep = ":"))

ran_xpnsl_windows =
	seq(1,29) %>%
	  future_map(
	    ~genwin(ran_xpnsl, .x) %>%
	      mutate(CHR = .x,
	             BP = (WindowStart + WindowStop)/2,
	             SNP = paste(CHR, BP, sep = ":")) %>%
	      select(SNP,CHR, BP, everything())
	  ) %>%
	  reduce(bind_rows)

ran_xpnsl %>%
  mutate(window = map2_chr(.x = CHR, .y = BP,
                           ~ windowfinder(.x, .y, ran_xpnsl_windows)[1])) %>%
	write_tsv("200910_SIM.xpnsl.out.norm.windows.txt")
