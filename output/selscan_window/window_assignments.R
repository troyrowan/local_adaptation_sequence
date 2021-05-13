library(furrr)
library(GenWin)
library(tidyverse)


windowfinder = function(chromosome, basepair, windows){
	window_id =
		windows %>%
		filter(CHR == chromosome &
					WindowStart <= basepair &
					WindowStop >= basepair) %>%
		.$SNP
		return(window_id)
	}

# sim_nsl_windows =
#   read_tsv("200907_SIM.all.nsl.windows.txt",
#            col_type = cols(SNP = col_character()))
#
# sim_nsl =
#   read_tsv("200907_SIM.somesim.nsl.out.100bins.norm.gz",
#              col_names = c("SNP", "BP", "freq", "sl1", "sl0", "nsl", "norm_nsl", "sig")) %>%
#     mutate(CHR = as.numeric(str_split_fixed(SNP, ":", n = 4)[,1]),
#            SNP = paste(CHR, BP, sep = ":"))
#
#
# sim_nsl %>%
#   #sample_n(1000) %>%
#   mutate(window = map2_chr(.x = CHR, .y = BP,
#                            ~ windowfinder(.x, .y, sim_nsl_windows)[1])) %>%
# 	write_tsv("200907_SIM.all.nsl.out.100bins.norm.windows.txt")

###################################################################################

# pb_sim_nsl_windows =
#   read_tsv("200907_SIM.purebred.nsl.windows.txt",
#            col_type = cols(SNP = col_character()))
#
# pb_sim_nsl =
#   read_tsv("200907_SIM.purebred.nsl.out.100bins.norm.gz",
#              col_names = c("SNP", "BP", "freq", "sl1", "sl0", "nsl", "norm_nsl", "sig")) %>%
#     mutate(CHR = as.numeric(str_split_fixed(SNP, ":", n = 4)[,1]),
#            SNP = paste(CHR, BP, sep = ":"))
#
# pb_sim_nsl %>%
#   #sample_n(1000) %>%
#   mutate(window = map2_chr(.x = CHR, .y = BP,
#                            ~ windowfinder(.x, .y, pb_sim_nsl_windows)[1])) %>%
# 	write_tsv("200907_SIM.purebred.nsl.out.100bins.norm.windows.txt")


###################################################################################

ran_nsl_windows =
  read_tsv("200910_RAN.all.nsl.windows.txt",
           col_type = cols(SNP = col_character()))

ran_nsl =
  read_tsv("200910_RAN.all.nsl.out.100bins.norm.gz",
             col_names = c("SNP", "BP", "freq", "sl1", "sl0", "nsl", "norm_nsl", "sig")) %>%
    mutate(CHR = as.numeric(str_split_fixed(SNP, ":", n = 4)[,1]),
           SNP = paste(CHR, BP, sep = ":"))


ran_nsl %>%
  #sample_n(1000) %>%
  mutate(window = map2_chr(.x = CHR, .y = BP,
                           ~ windowfinder(.x, .y, ran_nsl_windows)[1])) %>%
	 write_tsv("200910_RAN.all.nsl.out.100bins.norm.windows.txt")
