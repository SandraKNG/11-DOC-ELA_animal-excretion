  #### fish nutrient and DOC excretion across a lake DOC gradient ####
  # This code was created by S. Klemet-N'Guessan in 2022 and 2023
  # R version 4.3.0
  
  ################################ FINAL FIGURES #################################
  # load libraries ----
  library(RColorBrewer)
  library(MetBrewer) # Met brewer color palette
  library(ggdist) # for stat_halfeye
  library(ggpubr) # for ggarrange function that allows to put all graphs on 1 page
  library(patchwork) # to arrange plots in one figure
  library(factoextra) # for pca

  # Set up prediction data ----
  pred_DOC_df <- excr %>% tidyr::expand(Trophic.position2,
                                     AmDOC = c(seq(min(AmDOC), max(AmDOC),
                                                   length = 100),
                                               rep(median(AmDOC), 100)))
  pred_lnRR_df <- excr.vol %>% tidyr::expand(AmDOC = c(seq(min(AmDOC), max(AmDOC),
                                                           length = 100),
                                                       rep(median(AmDOC), 100)))
  
  # Prediction models ----
  gamNDOC.pred <- fitted_values(gamNDOC, pred_DOC_df)
  gamPDOC.pred <- fitted_values(gamPDOC, pred_DOC_df)
  gamNPDOC.pred <- fitted_values(gamNPDOC, pred_DOC_df)
  gamlnRRN.pred <- fitted_values(gamlnRRN, pred_lnRR_df)
  gamlnRRP.pred <- fitted_values(gamlnRRP, pred_lnRR_df)
  
  # Plot ----
  # ..set up plotting parameters ----
  Trophic.labels <- c('carnivore', 'omnivore')
  Species.labels <- c('fathead minnow', 'pearl dace', 
                      'white sucker', 'yellow perch')
  DOC.labels <- c('low', 'medium', 'high')
  DOM.labels <- c('terrestrial \n humic-like', 'SR (molecular size)', 'β:α (freshness)', 'HIX (humic)', 
                  'ubiquitous \n humic-like','soil \n fulvic-like', 'protein-like', 
                  'micrboial \n humic-like','FI (source)', expression(SUVA[254]~(aromaticity)), 'terrestrial \n humic-like')
  lnRR.labels <- c('DOC', 'terrestrial \n humic-like DOM', 'soil \n fulvic-like DOM',
                   'microbial \n humic-like DOM','protein-like DOM')
  point.size = 1.5
  point.size2 = 1
  point.alpha = .4
  line.size = 0.5
  text.size = 3
  
  # create plotting functions
  plot_gam <- function(df, y) {
    ggplot(df, aes(x = AmDOC, y = fitted)) +  
      geom_point(data = excr, aes(x = AmDOC, y = y,
                                  colour = Trophic.position2), 
                 size = point.size, alpha = point.alpha) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = Trophic.position2),
                  alpha = .2) +
      geom_line(aes(colour = Trophic.position2), linewidth = line.size) +
      theme_classic(base_size = 10) +
      guides(colour = guide_legend(
        override.aes = list(size = point.size + .5)
        )) +
      scale_color_manual(name = 'Trophic position',
                         labels = Trophic.labels,
                         values = c("#dd5129",  "#43b284")) +
      scale_fill_manual(name = 'Trophic position',
                         labels = Trophic.labels,
                         values = c("#dd5129","#43b284")) +
      theme(legend.text = element_text(size = 9))
  }
  
  plot_boxplot <- function(x) {
    ggplot(excr.aov, 
           aes(x = DOC.level, y = log10(x), 
               fill = DOC.level, colour = DOC.level)) +
      stat_halfeye(adjust = .5, width = .6, .width = 0, justification = -.4,
                   alpha = .3, point_colour = NA, aes(fill = DOC.level)) +
      geom_boxplot(width = .25, size = .5, outlier.shape = NA, alpha = .2) +
      geom_point(size = point.size2, alpha = point.alpha, 
                 position = position_jitter(seed = 1, width = .1)) +
      scale_x_discrete(labels = DOC.labels) +
      theme_classic(base_size = 10) +
      scale_fill_manual(values = met.brewer('Greek', 3, direction = -1)) +
      scale_colour_manual(values = met.brewer('Greek', 3, direction = -1)) 
  }
  
  plot_nmds <- function(df){
    ggplot(df, aes(x = NMDS1, y = NMDS2, fill = Trophic.position2)) +
      geom_point(aes(colour = Trophic.position2, shape = Trophic.position2), 
                 size = point.size) +
      theme_classic(base_size = 10)
  }
  
  # Figure 1 ----
  # # N excretion
  NexcrDOC.p <- plot_gam(gamNDOC.pred, excr$massnorm.N.excr) +
    labs(x = '',
         y = expression(atop('Mass-normalized',
                             paste(N~excretion~(μg~N~g^-1~h^-1))))) +
    theme(axis.text.x = element_blank()) +
    scale_x_continuous(breaks = c(3, 5, 7, 9, 11))
  NexcrDOC.p
  
  # P excretion
  PexcrDOC.p <- plot_gam(gamPDOC.pred, excr$massnorm.P.excr) +
    labs(x = '',
         y = expression(atop('Mass-normalized',
                             paste(P~excretion~(μg~P~g^-1~h^-1))))) +
    theme(axis.text.x = element_blank()) +
    scale_x_continuous(breaks = c(3, 5, 7, 9, 11))
  PexcrDOC.p
  
  # N:P excretion
  NPexcrDOC.p <- plot_gam(gamNPDOC.pred, log10(excr$massnorm.NP.excr)) +
    labs(x = expression(DOC~(mg~C~L^-1)),
         y = expression(atop(Log[10]~'mass-normalized', 
                             paste('N:P excretion (molar)')))) +
    geom_hline(yintercept = log10(16), linetype = 'dashed', linewidth = line.size) +
    annotate("text", x = 8.7, y = 1.5, label = 'Redfield ratio', size = text.size) +
    scale_x_continuous(breaks = c(3, 5, 7, 9, 11))
  NPexcrDOC.p
  
  
  # combine plots ----
  ggarrange(NexcrDOC.p, 
            PexcrDOC.p, 
            NPexcrDOC.p, 
            label.x = 0.2, label.y = 1,
            nrow = 3, ncol = 1, align = 'hv', heights = 1,
            common.legend = T,
            labels = c('(a)', '(b)', '(c)'),
            font.label = list(size = 10))
  
  ggsave('tables_figures/final_tables_figures/Fig1.tiff', 
         width = 10, height = 16.2,
         units = 'cm', dpi = 600, compression = 'lzw')
  
  # Figure 2 ----
  # ...DOC excretion vs. DOC ----
  Cexcr.p <- plot_boxplot(excr.aov$massnorm.C.excr) +
    labs(x = '',
         y = expression(atop(Log[10]~'mass-normalized',
                             paste(DOC~excretion~(mg~C~g^-1~h^-1))))) +
    annotate("text", x = 2.5, y = 1, label = '*', size = text.size +.5) +
    geom_segment(x = 2, xend = 3, y = .9, yend = .9, linewidth = line.size) +
    geom_segment(x = 2, xend = 2, y = .9, yend = .8, linewidth = line.size) +
    geom_segment(x = 3, xend = 3, y = .9, yend = .8, linewidth = line.size) +
    theme(axis.text.x = element_blank())
  Cexcr.p
  
  # ...DOC:N excretion vs. DOC ----
  CNexcr.p <- plot_boxplot(excr.aov$massnorm.CN.excr) +
    labs(x = '',
         y = expression(atop(Log[10]~'mass-normalized',
                             paste(DOC:N~excretion~(molar))))) +
    theme(axis.text.x = element_blank()) +
    annotate("text", x = 2.5, y = -1.15, label = '**', size = text.size + .5) +
    geom_segment(x = 2, xend = 3, y = -1.25, yend = -1.25, linewidth = line.size) +
    geom_segment(x = 2, xend = 2, y = -1.25, yend = -1.35, linewidth = line.size) +
    geom_segment(x = 3, xend = 3, y = -1.25, yend = -1.35, linewidth = line.size)
  CNexcr.p
  
  # ...DOC:P excretion vs. DOC ----
  CPexcr.p <- plot_boxplot(excr.aov$massnorm.CP.excr) +
    labs(x = 'DOC',
         y = expression(atop(Log[10]~'mass-normalized',
                             paste(DOC:P~excretion~(molar))))) +
    theme(axis.title.x = element_text(vjust = -.4, hjust = 0.46))
  CPexcr.p
  
  # combine plots ----
  ggarrange(Cexcr.p, CNexcr.p, CPexcr.p,
            labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 10), label.x = 0.17, label.y = 1,
            nrow = 3, align = "v", legend = 'none')
  
  ggsave('tables_figures/final_tables_figures/Fig2.tiff', 
         width = 10, height = 16.2, 
         units = 'cm', dpi = 600, compression = 'lzw')
  
  # Figure 3 ----
  DOMexcr.p <- ggplot(excr.DOM, aes(x = type, y = massnorm.excr)) +
    geom_boxplot(outlier.shape = NA, linewidth = line.size) +
    geom_jitter(aes(colour = DOC.level, fill = DOC.level, shape = Trophic.position2), shape = 21, 
                size = point.size, stroke = .25, alpha = .6, width = .25) +
    labs(x = 'DOM optical parameters',
         y = expression(Mass-normalized~excretion~rates~(g^-1~h^-1))) +
    geom_hline(aes(yintercept = 0), linetype = 'dashed') +
    guides(colour = guide_legend(override.aes = list(size = point.size + .5))) +
    scale_colour_manual(name = 'DOC',
                      labels = DOC.labels,
                      values = met.brewer('Greek', 3, direction = -1)) +
    scale_fill_manual(name = 'DOC',
                       labels = DOC.labels,
                       values = met.brewer('Greek', 3, direction = -1)) +
    scale_x_discrete(labels = DOM.labels) +
    theme_bw(base_size = 10) +
    theme(axis.title.x = element_text(vjust = -.4),
          axis.title.y = element_text(vjust = .7),
          axis.text.x = element_text(angle = 45, hjust = 1.1))
  DOMexcr.p
  
  ggsave('tables_figures/final_tables_figures/Fig3.tiff', 
         width = 7, height = 5, 
         units = 'in', dpi = 600, compression = 'lzw')
  
  # Figure 4 ----
  # ...low DOC ----
  nmds.l.p <- plot_nmds(nmds.l.scores) +
    stat_ellipse(level = .95, aes(colour = Trophic.position2)) +
    scale_shape_manual(values = c(8, 16, 15, 17)) +
    scale_color_manual(values = c("#f5c34d", "#dd5129", "#43b284")) +
    scale_fill_manual(values = c("#f5c34d", "#dd5129", "#43b284")) +
    annotate("text", x = .65, y = 1.15,
             label = 'low DOC', colour = "#f5c34d", size = text.size) +
    annotate("text", x = -1.4, y = .45,
               label = 'omnivore', colour = "#43b284", size = text.size) +
    annotate("text", x = .4, y = .3,
               label = 'carnivore', colour = "#dd5129", size = text.size) +
    xlab('')
  nmds.l.p
  
  # ...medium DOC ----
  nmds.m.p <- plot_nmds(nmds.m.scores) +
    scale_shape_manual(values = c(8, 16)) +
    scale_color_manual(values = c("#8d1c06", "#dd5129")) +
    annotate("text", x = -1.6, y = .35, 
             label = 'medium DOC', colour = "#8d1c06", size = text.size) +
    annotate("text", x = .9, y = .1,
             label = 'carnivore', colour = "#dd5129", size = text.size) +
    xlab('')
  nmds.m.p

  # ...high DOC ----
  nmds.h.p <- plot_nmds(nmds.h.scores) +
    scale_shape_manual(values = c(8, 16)) +
    scale_color_manual(values = c("#3c0d03", "#dd5129")) +
    annotate("text", x = 2, y = .55, 
             label = 'high DOC', colour = "#3c0d03", size = text.size) +
    annotate("text", x = -.2, y = -.6,
             label = 'carnivore', colour = "#dd5129", size = text.size)
  nmds.h.p
  
  # combine plots ----
  ggarrange(nmds.l.p, nmds.m.p, nmds.h.p,
            labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 10), label.x = 0.13, label.y = 0.98,
            nrow = 3, align = "h", legend = 'none')
  ggsave('tables_figures/final_tables_figures/Fig4.tiff', 
         width = 10, height = 16, 
         units = 'cm', dpi = 600, compression = 'lzw')
  
  # Figure 5 ----
  x <- 7.3
  lnRRN.p <- ggplot(gamlnRRN.pred, aes(x = AmDOC, y = fitted)) +
    geom_point(data = excr.vol, aes(x = AmDOC, y = lnRR.N), 
               size = point.size, alpha = point.alpha) +
    geom_ribbon(aes(ymin = lower, ymax = upper), 
                alpha = .2) +
    geom_line(linewidth = line.size) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    labs(x = expression(DOC~(mg~C~L^-1)),
         y = 'lnRR N') +
    theme_classic(base_size = 10) +
    #theme(axis.text.x = element_blank()) +
    annotate("text", x = x, y = 5,
             label = 'volumetric excretion above ambient concentration',
             size = text.size) +
    annotate("text", x = x, y = -5,
             label = 'volumetric excretion below ambient concentration',
             size = text.size)
  lnRRN.p
  
  lnRRP.p <- ggplot(excr.vol, aes(x = AmDOC, y = lnRR.P)) +
    geom_point(size = point.size, alpha = point.alpha) +
    geom_hline(yintercept = lnRR.P.av, linetype = 'longdash') +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    labs(x = expression(DOC~(mg~C~L^-1)),
         y = 'lnRR P') +
    theme_classic(base_size = 10) +
    annotate("text", x = x, y = 7,
             label = 'volumetric excretion above ambient concentration',
             size = text.size) +
    annotate("text", x = x, y = -1,
             label = 'volumetric excretion below ambient concentration',
             size = text.size)
  lnRRP.p
  
  lnRRDOM.p <- ggplot(excr.vol.lg, 
                      aes(x = DOC.level, y = variable)) +
    labs(x = 'DOC',
         y = '') +
    scale_x_discrete(labels = DOC.labels) +
    scale_y_discrete(labels = lnRR.labels) +
    geom_tile(aes(fill = value), colour = "white") +
    theme_bw(base_size = 10) +
    #theme(legend.position = 'bottom') +
    scale_fill_gradient2(name = 'lnRR', midpoint = 0, 
                         mid = "#eee8d5", high = "#dc322f", low = "#268bd2") +
    theme(plot.title = element_text(size = 10, face = 'bold')) +
    ggtitle("(c)")
  lnRRDOM.p
  
  
  lnRRNP.p <- ggarrange(
    lnRRN.p,
    lnRRP.p,
    labels = c("(a)", "(b)"),
    ncol = 2,
    font.label = list(size = 10),
    label.x = 0.14,
    label.y = 1.01,
    align = 'v'
  )
  lnRRNP.p
  
  ggarrange(lnRRNP.p, lnRRDOM.p, #lnRRP.p,
            nrow = 2, ncol = 1,
            font.label = list(size = 10), label.x = 0.14, label.y = 1.05, 
            align = 'v', heights = c(1.2, 1))
  
  #(lnRRN.p/lnRRP.p)|lnRRDOM.p + plot_annotation(tag_levels = "a")
  
  ggsave('tables_figures/final_tables_figures/Fig5.tiff', 
         width = 10, height = 9, bg = 'white',
         units = 'cm', dpi = 600, compression = 'lzw', scale = 1.7)
  
  # Figure S1 ----
  pca.all.p <- fviz_pca_biplot(pca.all,
                               repel = T,
                               label = 'var',
                               title = '',
                               ggtheme = theme_bw()) +
    xlim(-5,5) +
    ylim(-3,3) +
    labs(x = 'PC1 (55.3%)',
         y = 'PC2 (21.4%)')
  pca.all.p
  
  ggsave('tables_figures/final_tables_figures/FigS1.tiff', 
         width = 11, height = 11, 
         units = 'cm', dpi = 600, compression = 'lzw')
  
  # Figure S2 ----
  # plot PCA
  pca.DOM.p <- fviz_pca_biplot(pca.DOM,
                           repel = T,
                           label = 'var',
                           col.ind = excr.pca$AmDOC,
                           title = '',
                           ggtheme = theme_bw()) +
    xlim(-5,5) +
    ylim(-3,3) +
    labs(x = 'PC1 (75.2%)',
         y = 'PC2 (13%)') +
    scale_colour_gradient2(name = "DOC \n(mg C/L)", high = "#3c0d03")
  pca.DOM.p
  
  ggsave('tables_figures/final_tables_figures/FigS2.tiff', 
         width = 12, height = 11, 
         units = 'cm', dpi = 600, compression = 'lzw')
  
  # Figure S3 ----
  ggplot(excr.pca, aes(AmDOC, Epi.chla)) +
    geom_point(aes(colour = Site.name), size = 4) +
    #geom_smooth(method = lm, colour = 'black') +
    xlab('DOC (mg C/L)') +
    theme_classic(base_size = 20) +
    annotate('text', x = 4.5, y = 7, 
             size = 5, label = 'cor = 0.9, p < 0.001')
  
  # export final tables ----
  write_csv(excr, "output/excr_final.csv")
  write_csv(excr.vol, "output/excr_vol_final.csv")
  write_csv(excr.ss, "output/excr_summary.csv")
  write_csv(excrtp.ss, "output/excr_summary_site_tp.csv")
  write_csv(excrDOM.ss, "output/DOMexcr_summary.csv")

  