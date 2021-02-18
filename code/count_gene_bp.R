library(dplyr)
library(tidyr)
library(purrr)
library(readr)


genes =
  read_tsv("genes.txt") %>%
  filter()

genes %>%
  select(chr, start_pos, end_pos) %>%
  mutate(bp = map2(.x = start_pos, .y = end_pos, .f = ~seq(.x,.y))) %>%
  unnest(bp) %>%
  select(chr, bp) %>%
  unique() %>%
  nrow()
