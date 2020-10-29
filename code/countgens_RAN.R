library(pedigree)
library(tidyverse)

ran_ped =
  read_csv("data/200910_RAN/All_Animals_SireDam.csv") %>%
  select(id = anm_key, sire = sire_key, dam = dam_key)

gens = countGen(ran_ped)


bind_cols(ran_ped, gens) %>%
  write_csv("RAN_ped_countGens.csv")
