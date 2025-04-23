  #### Fish supply distinct nutrients and dissolved organic matter composition ####
  ### relative to ambient concentrations in northern lakes ####
  
  # This code was created by S. Klemet-N'Guessan in 2022-2024
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
  pred_df <- excr %>% tidyr::expand(Trophic.position2, 
                                     AmDOC = c(seq(min(AmDOC), max(AmDOC),
                                                   length = 100),
                                               rep(median(AmDOC), 100)))
  
  # Prediction models ----
  gamNDOC.pred <- fitted_values(gamNDOC, pred_df)
  gamPDOC.pred <- fitted_values(gamPDOC, pred_df)
  gamNPDOC.pred <- fitted_values(gamNPDOC, pred_df)
  
  # Plot ----
  # ..set up plotting parameters ----
  Trophic.labels <- c('carnivore', 'omnivore')
  Species.labels <- c('fathead minnow', 'northern pearl dace', 
                      'white sucker', 'yellow perch')
  DOC.labels <- c('low', 'medium', 'high')
  DOM.labels <- c('C3 (terrestrial \n humic-like)', 'C1 (ubiquitous \n humic-like)', 
                  'HIX (humic)', 'C4 (soil \n fulvic-like)', 'SR (molecular size)', 
                  'β:α (freshness)', 'C7 (protein-like)', 'C5 (microbial \n humic-like)',
                  'C2 (terrestrial \n humic-like)', expression(SUVA[254]~(aromaticity)), 
                  'FI (source)')
  point.size = 1.5
  point.size2 = 1
  point.alpha = .4
  line.size = 0.5
  text.size = 3
  text.size2 = 2.5
  x = 7.3
  
  # create plotting functions
  plot_gam <- function(df, x, y) {
    ggplot(df, aes(x = .data[[x]], y = fitted)) +
      geom_point(data = excr, aes(x = .data[[x]], y = y,
                                  colour = Trophic.position2),
                 size = point.size, alpha = point.alpha) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = Trophic.position2),
                  alpha = .2) +
      geom_line(linewidth = line.size, aes(colour = Trophic.position2)) +
      scale_x_continuous(n.breaks = 6) +
      theme_classic(base_size = 10) +
      guides(colour = guide_legend(
        override.aes = list(size = point.size + .5)
        )) +
      scale_color_manual(name = 'Trophic guild',
                         labels = Trophic.labels,
                         values = c("#dd5129",  "#43b284")) +
      scale_fill_manual(name = 'Trophic guild',
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
    ggplot(df, aes(x = NMDS1, y = NMDS2, fill = Trophic.position2, shape = Trophic.position2)) +
      geom_point(aes(colour = Trophic.position2), 
                 size = point.size) +
      theme_classic(base_size = 10)
  }
  
  plot_coeff <- function(x, df, y1, y2) {
    ggplot(df, aes(x = log10(Mass), y = log10(.data[[x]]))) +
      geom_point(size = 0.9, alpha = point.alpha) +
      geom_smooth(method = 'lm', color = 'black') +
      xlab('') +
      stat_regline_equation(label.y = y1, aes(label = after_stat(eq.label)), 
                            size = text.size2) +
      stat_regline_equation(label.y = y2, aes(label = ..rr.label..), 
                            size = text.size2) 
      theme_classic(base_size = 8)
  }
  
  plot_sp <- function(x, df){
    ggplot(df,
           aes(x = AmDOC, y = .data[[x]])) +
      geom_point(aes(colour = Species.code),
                 size = point.size + .2,
                 alpha = point.alpha) +
      theme_classic(base_size = 10) +
      scale_x_continuous(n.breaks = 6) +
      scale_color_brewer(palette = 'Dark2',
                         name = 'Species',
                         labels = Species.labels)
  }
  
  plot_sp_bp <- function(x, df){
    
  }
   
  
  # Figure 1 ----
  # ...DOC ----
  # N excretion
  NexcrDOC.p <- plot_gam(gamNDOC.pred, 'AmDOC', excr$massnorm.N.excr) +
    labs(x = '',
         y = expression(atop('Mass-normalized',
                             paste(N~excretion~(μg~N~g^-1~h^-1)))))
  NexcrDOC.p
  
  # P excretion
  PexcrDOC.p <- plot_gam(gamPDOC.pred %>% filter(Trophic.position2 == 'C'), 'AmDOC',
                         excr$massnorm.P.excr) +
    labs(x = '',
         y = expression(atop('Mass-normalized',
                             paste(P~excretion~(μg~P~g^-1~h^-1))))) +
    geom_hline(data = excrtp.ss %>% filter(Variable == 'massnorm.P.excr' &
                                           .group == 'Trophic.position2=O'),
               aes(yintercept = Mean), linetype = 'dashed',
               linewidth = line.size, colour = "#43b284")
  PexcrDOC.p
  
  # N:P excretion vs DOC
  NPexcrDOC.p <- plot_gam(gamNPDOC.pred, 'AmDOC', log10(excr$massnorm.NP.excr)) +
    labs(x = expression(DOC~(mg~C~L^-1)),
         y = expression(atop(Log[10]~'mass-normalized', 
                             paste('N:P excretion (molar)')))) +
    geom_hline(yintercept = log10(16), linewidth = line.size) +
    annotate("text", x = 4.5, y = 1.5, label = 'Redfield ratio', size = text.size) 
  NPexcrDOC.p
  
  #...Zmean ----
  # N excretion
  NexcrZmean.p <- plot_gam(gamNZmean.pred, 'Zmean', excr$massnorm.N.excr) +
    labs(x = '',
         y = '') 
  NexcrZmean.p
  
  # P excretion
  PexcrZmean.p <- plot_gam(gamPZmean.pred, 'Zmean', #%>% filter(Trophic.position2 == 'C'), 
                           excr$massnorm.P.excr) +
    labs(x = '',
         y = '') #+
  # geom_hline(data = excrtp.ss %>% filter(Variable == 'massnorm.P.excr' &
  #                                        .group == 'Trophic.position2=O'),
  #            aes(yintercept = Mean), linetype = 'dashed', 
  #            linewidth = line.size, colour = "#43b284")
  PexcrZmean.p
  
  # N:P excretion vs Zmean
  NPexcrZmean.p <- plot_gam(gamNPZmean.pred, 'Zmean', log10(excr$massnorm.NP.excr)) +
    labs(x = 'Mean depth (m)',
         y = '') +
    geom_hline(yintercept = log10(16), linewidth = line.size) +
    annotate("text", x = 4.5, y = 1.5, label = 'Redfield ratio', size = text.size) 
  NPexcrZmean.p
  
  # combine plots ----
  ggarrange(NexcrDOC.p, #NexcrZmean.p, 
            PexcrDOC.p, #PexcrZmean.p,
            NPexcrDOC.p, #NPexcrZmean.p,
            label.x = 0.2, label.y = 1,
            nrow = 3, ncol = 1, align = 'hv',
            common.legend = T,
            labels = c('(a)', '(b)', '(c)'),#, '(d)', '(e)', '(f)'),
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
    geom_jitter(aes(colour = DOC.level, fill = DOC.level, shape = Trophic.position2),  
                size = point.size, stroke = .25, alpha = .6, width = .25) +
    labs(x = 'DOM optical parameters',
         y = expression(Mass-normalized~excretion~rates~(g^-1~h^-1))) +
    guides(colour = guide_legend(override.aes = list(size = point.size + .5))) +
    scale_colour_manual(name = 'DOC',
                      labels = DOC.labels,
                      values = met.brewer('Greek', 3, direction = -1)) +
    scale_fill_manual(name = 'DOC',
                       labels = DOC.labels,
                       values = met.brewer('Greek', 3, direction = -1)) +
    scale_shape_manual(name = 'Trophic guild',
                      labels = Trophic.labels,
                      values = c(16, 17)) +
    scale_x_discrete(labels = DOM.labels) +
    theme_bw(base_size = 10) +
    theme(axis.title.x = element_text(vjust = -.4),
          axis.title.y = element_text(vjust = .7),
          axis.text.x = element_text(angle = 45, hjust = 1.1))
  DOMexcr.p
  
  ggsave('tables_figures/final_tables_figures/Fig3_2.tiff', 
         width = 7, height = 5, 
         units = 'in', dpi = 600, compression = 'lzw', scale = .9)
  
  ggsave('tables_figures/final_tables_figures/Fig3_2.png', 
         width = 7, height = 5, 
         units = 'in', dpi = 600,  scale = .9)
  
  # Figure 4 ----
  # ...low DOC ----
  nmds.l.p <- plot_nmds(nmds.l.scores) +
    scale_shape_manual(values = c(8, 16, 16)) +
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
            nrow = 3, align = "v", legend = 'none')
  ggsave('tables_figures/final_tables_figures/Fig4.tiff', 
         width = 10, height = 16, 
         units = 'cm', dpi = 600, compression = 'lzw')
  
  ggsave('tables_figures/final_tables_figures/Fig4_2.png', 
         width = 10, height = 18, 
         units = 'cm', dpi = 600)
  
  
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
    scale_colour_gradient2(name = expression(atop(DOC,
                                                  paste((mg~C~L^-1)))), high = "#3c0d03")
  pca.DOM.p
  
  ggsave('tables_figures/final_tables_figures/FigS2.tiff', 
         width = 12, height = 11, 
         units = 'cm', dpi = 600, compression = 'lzw')
  
  # Figure S3 ----
  coeffN.p <- plot_coeff("N.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~N~excretion, paste((μg~N~ind^-1~h^-1))))) +
    stat_regline_equation(label.y = 3, aes(label = ..rr.label..), 
                          size = text.size2)
  
  coeffP.p <- plot_coeff("P.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~P~excretion, paste((μg~P~ind^-1~h^-1))))) +
    stat_regline_equation(label.y = 2.7, aes(label = ..rr.label..), 
                          size = text.size2) 
  
  coeffC.p <- plot_coeff("C.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~DOC~excretion, paste((mg~C~ind^-1~h^-1))))) +
    stat_regline_equation(label.y = 1, aes(label = ..rr.label..), 
                          size = text.size2) 
  
  coeffSUVA.p <- plot_coeff("SUVA.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~SUVA[254]~excretion, paste((L~mg~C~m^-1~ind^-1~h^-1))))) +
    stat_regline_equation(label.y = 0.1, aes(label = ..rr.label..), 
                          size = text.size2) 
  
  coeffSR.p <- plot_coeff("SR.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~SR~excretion, paste((ind^-1~h^-1))))) +
    stat_regline_equation(label.y = -0.6, aes(label = ..rr.label..), 
                          size = text.size2) 
  
  coeffBA.p <- plot_coeff("BA.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~β:α~excretion, paste((ind^-1~h^-1))))) +
    stat_regline_equation(label.y = -0.8, aes(label = ..rr.label..), 
                          size = text.size2) 
  
  coeffFI.p <- plot_coeff("FI.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~FI~excretion, paste((ind^-1~h^-1))))) +
    stat_regline_equation(label.y = -0.1, aes(label = ..rr.label..), 
                          size = text.size2) 
  
  coeffHIX.p <- plot_coeff("HIX.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~HIX~excretion, paste((ind^-1~h^-1))))) +
    stat_regline_equation(label.y = -0.6, aes(label = ..rr.label..), 
                          size = text.size2) 
  
  coeffC1.p <- plot_coeff("C1.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~C1~excretion, paste((RU~ind^-1~h^-1))))) +
    stat_regline_equation(label.y = -0.9, aes(label = ..rr.label..), 
                          size = text.size2) 
  
  coeffC2.p <- plot_coeff("C2.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~C2~excretion, paste((RU~ind^-1~h^-1))))) +
    stat_regline_equation(label.y = -0.1, aes(label = ..rr.label..), 
                          size = text.size2) 
  
  coeffC3.p <- plot_coeff("C3.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~C3~excretion, paste((RU~ind^-1~h^-1))))) +
    stat_regline_equation(label.y = -1.4, aes(label = ..rr.label..), 
                          size = text.size2) 
  
  coeffC4.p <- plot_coeff("C4.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~C4~excretion, paste((RU~ind^-1~h^-1))))) +
    stat_regline_equation(label.y = -1, aes(label = ..rr.label..), 
                          size = text.size2) 
  
  coeffC5.p <- plot_coeff("C5.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~C5~excretion, paste((RU~ind^-1~h^-1))))) +
    stat_regline_equation(label.y = -0.3, aes(label = ..rr.label..), 
                          size = text.size2) 
  
  coeffC7.p <- plot_coeff("C7.excretion.rate", excr.var) +
    ylab(expression(atop(Log[10]~C7~excretion, paste((RU~ind^-1~h^-1))))) +
    stat_regline_equation(label.y = -0.2, aes(label = ..rr.label..), 
                          size = text.size2) 
  
  
  # combine plots ----
  figS2 <- ggarrange(coeffN.p, coeffP.p, coeffC.p, coeffSUVA.p, coeffSR.p, 
            coeffBA.p, coeffFI.p, coeffHIX.p, coeffC1.p, coeffC2.p, 
            coeffC3.p, coeffC4.p, coeffC5.p, coeffC7.p,
            labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)", 
                      "(h)", "(i)", "(j)", "(k)", "(l)", "(m)", "(n)"),
            font.label = list(size = 10), label.x = 0.13, label.y = 1.02, 
            align = "hv", legend = 'none')
  annotate_figure(figS2, 
                  bottom = text_grob(expression(Log[10]~mass~(g)), size = 10, y = 1))
  ggsave('tables_figures/final_tables_figures/FigS3.tiff', 
         width = 18, height = 19, 
         units = 'cm', dpi = 600, compression = 'lzw', bg = 'white')
  
  # Figure S4 ----
  pca.all.p <- fviz_pca_biplot(pca.all,
                               repel = T,
                               label = 'var',
                               title = '',
                               ggtheme = theme_bw()) +
    xlim(-5,5) +
    ylim(-3,3) +
    labs(x = 'PC1 (58.4%)',
         y = 'PC2 (18.7%)')
  pca.all.p
  
  pca.result
  
  ggsave('tables_figures/final_tables_figures/FigS2.tiff', 
         width = 11, height = 11, 
         units = 'cm', dpi = 600, compression = 'lzw')
  
  
  # Figure S6 ----
  # Individual excretion rates across streams by species
  Nexcrsp.p <- plot_sp('massnorm.N.excr', excr) +
    labs(x = '',
         y = expression(atop('Mass-normalized',
                             paste(N~excretion~(μg~N/g/h)))))
  Nexcrsp.p
  
  Pexcrsp.p <- plot_sp('massnorm.P.excr', excr) +
    labs(x = '',
         y = expression(atop('Mass-normalized',
                             paste(P~excretion~(μg~P/g/h))))) 
  Pexcrsp.p 
  
  NPexcrsp.p <- plot_sp('log10.massnorm.NP.excr', excr) +
    labs(x = expression(DOC~(mg~C~L^-1)),
         y = expression(atop(Log[10]~'mass-normalized', 
                             paste('N:P excretion (molar)')))) 
  NPexcrsp.p 
  
  # boxplots of species variability
  ggplot(excr, aes(x = Species.code, y = massnorm.N.excr)) +
    geom_point(aes(colour = Species.code)) +
    geom_boxplot() +
    scale_color_brewer(palette = 'Dark2',
                       name = 'Species',
                       labels = Species.labels)
  
  # combine plots ----
  ggarrange(Nexcrsp.p, Pexcrsp.p, NPexcrsp.p,
            labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 10), label.x = 0.17, label.y = 1,
            nrow = 3, align = "v", common.legend = T, legend = 'right')
  
  ggsave('tables_figures/final_tables_figures/FigS6.tiff', 
         width = 6, height = 6, 
         units = 'in', dpi = 600, compression = 'lzw')
  
  # Figure S7 ----
  # NMDS of all ambient DOC + fish excr
  nmds.all.p <- ggplot(nmds.all.scores, aes(x = NMDS1, y = NMDS2, 
                                            fill = Source, shape = Source)) +
    geom_point(aes(colour = Source), 
               size = point.size) +
    theme_classic(base_size = 10) +
    scale_shape_manual(values = c(8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 16, 16, 16, 16)) +
    scale_color_manual(values = c("grey", "#3c0d03", "#f5c34d", "#8d1c06",
                                  "grey", "grey", "grey",
                                  "grey", "grey", "grey", "grey",
                                  'darkgoldenrod', 'salmon', 'bisque3', 'olivedrab')) +
    scale_fill_manual(values = c("grey", "#3c0d03", "#f5c34d", "#8d1c06",
                                 "grey", "grey", "grey",
                                 "grey", "grey", "grey", "grey",
                                 'darkgoldenrod', 'salmon', 'bisque3', 'olivedrab')) +
    theme(legend.position = "none") +
    annotate("text", x = 1, y = .7,
             label = 'low DOC', colour = "#f5c34d", size = text.size) +
    annotate("text", x = 1.77, y = .7,
             label = 'high DOC', colour = "#3c0d03", size = text.size) +
    annotate("text", x = 1.6, y = 0,
             label = 'medium DOC', colour = "#8d1c06", size = text.size) +
  annotate("text", x = -.5, y = .9,
           label = 'fathead minnow', colour = 'darkgoldenrod', size = text.size) +
    annotate("text", x = .5, y = -.75,
             label = 'northern pearl dace', colour = 'salmon', size = text.size) +
    annotate("text", x = .5, y = 0.2,
             label = 'white sucker', colour = 'bisque3', size = text.size) +
      annotate("text", x = -1.5, y = -1,
             label = 'yellow perch', colour = 'olivedrab', size = text.size) 
  nmds.all.p
  
  ggsave('tables_figures/final_tables_figures/FigS7_sp.tiff', 
         width = 12, height = 8, 
         units = 'cm', dpi = 600, compression = 'lzw')
  
  # export final tables ----
  write_csv(excr, "output/excr_final.csv")
  write_csv(excr.vol, "output/excr_vol_surf_final2.csv")
  write_csv(excr.ss, "output/excr_summary.csv")
  write_csv(excrtp.ss, "output/excr_summary_site_tp.csv")
  write_csv(excr.DOM.ss, "output/DOMexcr_summary.csv")
  write_csv(excr.pca, "output/AmDOM_summary.csv")
  write_csv(fish.biomass, "output/fish_biomass.csv")
  write_csv(biomass.av, "output/fish_biomass_av.csv")

  