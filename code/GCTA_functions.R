read_blp = function(blp_path){
  blp_results =
    read_tsv(blp_path,
             col_names = c("FID", "international_id", "IV1", "BV", "IV2", "Residual")) %>%
    select(international_id, BV, Residual)
  return(blp_results)
}
library(qvalue)
library(tidyverse)

read_gwas = function(gwas_path){
  gwas_results =
    read_tsv(gwas_path) %>%
    filter(p < 0.01) %>%
    mutate(chrbp = paste(Chr, bp, sep = ":"),
           maf = case_when(Freq > 0.5 ~ 1-Freq,
                           TRUE ~ Freq)) %>%
    dplyr::rename(CHR = Chr, BP = bp)
  return(gwas_results)
}

make_qq <- function(dd, x) {
  dd<-dd[order(dd[[x]]), ]
  dd$qq <- qnorm(ppoints(nrow(dd)))
  dd
}

ggmanhattan2 = function(inputfile,
                        prune = 1.0,
                        value = p,
                        alpha = 0.5,
                        pcol = "p",
                        sig_threshold_p = 1e-5,
                        sig_threshold_q = 0.1,
                        pointsize = 1.0,
                        colors = c("grey10", "grey55"),
                        sigsnps = NULL,
                        sigsnps2 = NULL,
                        sigsnps3 = NULL,
                        sigsnpcolor = "springgreen3",
                        sigsnpcolor2 = "goldenrod1",
                        sigsnpcolor3 = "royalblue4"){
  require(qvalue)
  require(dplyr)
  require(stringr)
  require(ggplot2)
  require(cowplot)
  require(viridis)
  #Allows p or q value to be interpreted as column name in ggplot call
  v = enexpr(value)
  gwas = inputfile
  #print(v)
  # gwas = inputfile %>% #reads in input file
  #   mutate(q = qvalue(p)$qvalues) %>% #transforms p-values to q-values
  #   #ifelse(v == "p", filter(p<prune), filter(q<prune))
  #   filter(., if (v == "p") p < prune else q < prune)

  don = gwas %>%
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    left_join(gwas, ., by = "CHR") %>%
    arrange(CHR, BP) %>%
    mutate(BPcum=BP+tot) %>%
    ungroup()

  axisdf = don %>%
    group_by(CHR) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  gwas_plot =
    don %>%
    ggplot(
      aes(x=BPcum, y=-log10(!!v))) +
    #geom_point(aes(color=as.factor(chr)), alpha=alpha, size=1.3) +
    geom_point(aes(color=as.factor(CHR)), alpha=alpha, size = pointsize) +
    scale_color_manual(values = rep(colors, 29)) +
    scale_x_continuous(label = axisdf$CHR[c(TRUE, FALSE, FALSE)], breaks= axisdf$center[c(TRUE, FALSE, FALSE)] ) +
    #scale_y_continuous(expand = c(0, 0.01) ) +
    labs(x = "",
         y = case_when(v == "p" ~ expression(paste(-log10, "(", italic('p'), ")")),
                       v == "q" ~ expression(paste(-log10, "(", italic('q'), ")"))))+
    theme_bw() +
    geom_point(data = subset(don, SNP %in% sigsnps), color=sigsnpcolor, size = pointsize, alpha = 1) +
    geom_point(data = subset(don, SNP %in% sigsnps2), color=sigsnpcolor2, size = pointsize, alpha = 1) +
    geom_point(data = subset(don, SNP %in% sigsnps3), color=sigsnpcolor3, size = pointsize, alpha = 1) +
    geom_hline(yintercept = case_when(v == "p" ~ -log10(sig_threshold_p),
                                      v == "q" ~ -log10(sig_threshold_q)),
               color = "red",
               size = 0.25) +
    theme(legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()
    )+
    theme(
      panel.background = element_rect(fill = "transparent") # bg of the panel
      , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
      , panel.grid.major = element_blank() # get rid of major grid
      , panel.grid.minor = element_blank() # get rid of minor grid
      , legend.background = element_rect(fill = "transparent") # get rid of legend bg
      , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )#+
  # theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
  #       axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0.5, face = "plain"),
  #       axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "bold"),
  #       axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = 2.5, face = "bold"))
  return(gwas_plot)
}



selscan_manhattans = function(inputfile,
                        col = ihs,
                        alpha = 0.5,
                        pointsize = 1.0,
                        colors = c("grey10", "grey55"),
                        ylab = "",
                        sigsnps = NULL,
                        sigsnps2 = NULL,
                        sigsnps3 = NULL,
                        sigsnpcolor = "springgreen3",
                        sigsnpcolor2 = "goldenrod1",
                        sigsnpcolor3 = "royalblue4"){
  require(qvalue)
  require(dplyr)
  require(stringr)
  require(ggplot2)
  require(cowplot)
  require(viridis)
  #Allows p or q value to be interpreted as column name in ggplot call
  v = enexpr(col)
  #print(v)
  gwas = inputfile

  don = gwas %>%
    group_by(CHR) %>%
    summarise(chr_len=max(BP)) %>%
    mutate(tot=cumsum(chr_len)-chr_len) %>%
    select(-chr_len) %>%
    left_join(gwas, ., by = "CHR") %>%
    arrange(CHR, BP) %>%
    mutate(BPcum=BP+tot) %>%
    ungroup()

  axisdf = don %>%
    group_by(CHR) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  gwas_plot =
    don %>%
    ggplot(
      aes(x=BPcum, y=!!v)) +
    #geom_point(aes(color=as.factor(chr)), alpha=alpha, size=1.3) +
    geom_point(aes(color=as.factor(CHR)), alpha=alpha, size = pointsize) +
    scale_color_manual(values = rep(colors, 29)) +
    scale_x_continuous(label = axisdf$CHR[c(TRUE, FALSE, FALSE)], breaks= axisdf$center[c(TRUE, FALSE, FALSE)] ) +
    #scale_y_continuous(expand = c(0, 0.01) ) +
    labs(x = "",
         y = "")+
    theme_bw() +
    geom_point(data = subset(don, SNP %in% sigsnps), color=sigsnpcolor, size = pointsize, alpha = 0.5) +
    geom_point(data = subset(don, SNP %in% sigsnps2), color=sigsnpcolor2, size = pointsize, alpha = 0.5) +
    geom_point(data = subset(don, SNP %in% sigsnps3), color=sigsnpcolor3, size = pointsize, alpha = 0.5) +
    theme(legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()
    )+
    theme(
      panel.background = element_rect(fill = "transparent") # bg of the panel
      , plot.background = element_rect(fill = "transparent", color = NA) # bg of the plot
      , panel.grid.major = element_blank() # get rid of major grid
      , panel.grid.minor = element_blank() # get rid of minor grid
      , legend.background = element_rect(fill = "transparent") # get rid of legend bg
      , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )#+
  # theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
  #       axis.text.y = element_text(color = "grey20", size = 20, angle = 0, hjust = 1, vjust = 0.5, face = "plain"),
  #       axis.title.x = element_text(color = "grey20", size = 20, angle = 0, hjust = .5, vjust = 0, face = "bold"),
  #       axis.title.y = element_text(color = "grey20", size = 20, angle = 90, hjust = .5, vjust = 2.5, face = "bold"))
  return(gwas_plot)
}


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

StatQq2 <- ggproto("StatQq", Stat,
                   default_aes = aes(y = after_stat(sample), x = after_stat(theoretical)),

                   required_aes = c("sample"),

                   compute_group = function(data, scales, quantiles = NULL,
                                            distribution = stats::qnorm, dparams = list(),
                                            na.rm = FALSE) {

                     sample <- sort(data$sample)
                     n <- length(sample)

                     # Compute theoretical quantiles
                     if (is.null(quantiles)) {
                       quantiles <- stats::ppoints(n)
                     } else if (length(quantiles) != n) {
                       abort("length of quantiles must match length of data")
                     }

                     theoretical <- do.call(distribution, c(list(p = quote(quantiles)), dparams))

                     res <- ggplot2:::new_data_frame(list(sample = sample,
                                                          theoretical = theoretical))

                     # NEW: append remaining columns from original data
                     # (e.g. if there were other aesthetic variables),
                     # instead of returning res directly
                     data.new <- subset(data[rank(data$sample), ],
                                        select = -c(sample, PANEL, group))
                     if(ncol(data.new) > 0) res <- cbind(res, data.new)

                     res
                   }
)

stat_qq2 <- function (mapping = NULL, data = NULL, geom = "point",
                      position = "identity", ..., distribution = stats::qnorm,
                      dparams = list(), na.rm = FALSE, show.legend = NA,
                      inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = StatQq2, geom = geom,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(distribution = distribution, dparams = dparams,
                      na.rm = na.rm, ...))
}


