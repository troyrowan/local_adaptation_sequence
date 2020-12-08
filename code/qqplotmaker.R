library(tidyverse)

library(cowplot)


ggqq <- function(pvector){
  pvector = pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
  pdf = data.frame(observed = -log10(sort(pvector,decreasing=FALSE)), expected = -log10(ppoints(length(pvector))))
  qqplotted = ggplot(pdf, aes(expected, observed))+
    geom_point()+
    geom_abline(intercept = 0,
                slope = 1,
                colour = "red")+
    labs(x = expression(paste("Expected ", -log10, "(", italic('p'), ")")),
         y = expression(paste("Observed ", -log10, "(", italic('p'), ")")))+
    theme_cowplot()
  return(qqplotted)
}

plot_grid(
  read_tsv("output/200910_RAN/seq_gwas/200910_RAN.full_age.pvalues.txt.gz",
           col_names = c("p"),
           col_types = cols(p = col_double())) %>%
    .$p %>%
    ggqq(),
  read_tsv("output/200910_RAN/seq_gwas/200910_RAN.young_age.pvalues.txt.gz",
           col_names = c("p"),
           col_types = cols(p = col_double())) %>%
    .$p %>%
    ggqq(),
  read_tsv("output/200907_SIM/seq_gwas/200907_SIM.full_age.pvalues.txt.gz",
           col_names = c("p"),
           col_types = cols(p = col_double())) %>%
    #filter(p < 0.1) %>%
    .$p %>%
    ggqq(),
  read_tsv("output/200907_SIM/seq_gwas/200907_SIM.pb_age.pvalues.txt.gz",
           col_names = c("p"),
           col_types = cols(p = col_double())) %>%
    .$p %>%
    ggqq(),
  nrow = 2,
  labels = "AUTO"
  )

ggsave("output/figures/QQPLOTS.png", height = 10, width = 10, units = "in")
