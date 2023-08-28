  #### fish nutrient and DOC excretion across a lake DOC gradient ####
  # This code was created by S. Klemet-N'Guessan in 2022 and 2023
  # R version 4.3.0
  
  ################################ FINAL FIGURES #################################
  # load libraries ----
  library(RColorBrewer)
  library(MetBrewer) # Met brewer color palette
  library(ggdist) # for stat_halfeye
  library(ggpubr) # for ggarrange function that allows to put all graphs on 1 page
  library(factoextra) # for pca

  
  # Set up prediction data ----
  # Prediction function for trophic position
  # make the prediction, add this and a column of standard errors to the 
  # prediction data.frame.
  
  pred_df <- function(x) {
    newdf <- with(excr, 
                  expand.grid(x = c(seq(min(x), max(x),
                                        length = 100), 
                                    rep(median(x), 100))))
    
    function(mod) {
      newdf_pred <- cbind(newdf, 
                          predict(mod, newdata = newdf,
                                  se.fit = TRUE,
                                  type = "link"))
      
      ilink <- inv_link(mod)
      newdf_pred <- newdf_pred %>%
        mutate(response = ilink(fit),
               upper = ilink(fit + (1.96 * se.fit)),
               lower = ilink(fit - (1.96 * se.fit)))
      
      return(newdf_pred)
    }
  }
  
  # Prediction function for vertebrate classification with  pred on log scale 
  pred_log_df <- function(x) {
    newdf <- with(excr, 
                  expand.grid(x = c(seq(min(x), max(x),
                                        length = 100), 
                                    rep(median(x), 100))))
    
    function(mod) {
      newdf_pred <- cbind(newdf, 
                          predict(mod, newdata = newdf,
                                  se.fit = TRUE,
                                  type = "response"))
      
      newdf_pred <- newdf_pred %>% 
        mutate(upper = fit + (1.96 * se.fit),
               lower = fit - (1.96 * se.fit))
      
      return(newdf_pred)
    }
  }
  
  # Storing prediction functions with different inputs
  pred.DOC <- pred_df(excr$AmDOC)
  pred.DOM <- pred_df(excr$PC1)
  pred_log.DOC <- pred_log_df(excr$AmDOC)
  pred_log.DOM <- pred_log_df(excr$PC1)
  
  # Prediction models ----
  gamNDOC.pred <- pred.DOC(gamNDOC)
  gamPDOC.pred <- pred.DOC(gamPDOC)
  gamNDOM.pred <- pred.DOM(gamNDOM)
  gamPDOM.pred <- pred.DOM(gamPDOM)
  gamNPDOC.pred <- pred_log.DOC(gamNPDOC)
  gamNPDOM.pred <- pred_log.DOM(gamNPDOM)
  
  # Plot ----
  # ..set up plotting parameters ----
  Trophic.labels <- c('invert/piscivore', 'invertivore', 'omnivore')
  Species.labels <- c('fathead minnow', 'pearl dace', 
                      'white sucker', 'yellow perch')
  DOC.labels <- c('low', 'medium', 'high')
  point.size = 1.5
  point.size2 = 1
  point.alpha = .4
  line.size = 0.5
  text.size = 3
  
  # create plotting functions
  plot_gam <- function(df) {
    ggplot(df, aes(x = x, y = response)) +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = .2) +
      geom_line(linewidth = line.size) +
      theme_classic(base_size = 10) +
      guides(colour = guide_legend(
        override.aes = list(size = point.size + .5)
        )) +
      scale_color_manual(name = 'Trophic position',
                         labels = Trophic.labels,
                         values = c("#dd5129", "steelblue4", "#43b284")) +
      theme(legend.text = element_text(size = 9))
  }
  
  plot_log_gam <- function(df) {
    ggplot(df, aes(x = x, y = fit)) +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = .2) +
      geom_line(linewidth = line.size) +
      theme_classic(base_size = 10) +
      guides(colour = guide_legend(override.aes = 
                                     list(size = point.size + .5))) +
      scale_color_manual(name = 'Trophic position',
                         labels = Trophic.labels,
                         values = c("#dd5129", "steelblue4", "#43b284")) 
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
    ggplot(df, aes(x = NMDS1, y = NMDS2, fill = Trophic.position)) +
      geom_point(aes(colour = Trophic.position, shape = Trophic.position), 
                 size = point.size) +
      theme_classic(base_size = 10)
  }
  
  # Figure 1 ----
  # ...DOC ----
  # N excretion
  NexcrDOC.p <- plot_gam(gamNDOC.pred) +
    geom_point(data = excr, aes(x = AmDOC, y = massnorm.N.excr,
                                   colour = Trophic.position), 
               size = point.size, alpha = point.alpha) +
    labs(x = '',
         y = expression(atop('Mass-normalized',
                             paste(N~excretion~(μg~N/g/h))))) +
    theme(axis.text.x = element_blank()) +
    scale_x_continuous(breaks = c(3, 5, 7, 9, 11))
  NexcrDOC.p
  
  # P excretion
  PexcrDOC.p <- ggplot(excr, aes(x = AmDOC, y = massnorm.P.excr,
                                           colour = Trophic.position)) +
    geom_point(size = point.size, alpha = point.alpha) +
    geom_hline(data = excr.ss %>%
                 filter(Variable == 'massnorm.P.excr'), aes(yintercept = Mean),
                linetype = 'dashed', linewidth = line.size) +
    labs(x = '',
         y = expression(atop('Mass-normalized',
                             paste(P~excretion~(μg~P/g/h))))) +
    theme_classic(base_size = 10) +
    scale_color_manual(name = 'Trophic position',
                       labels = Trophic.labels,
                       values = met.brewer('Egypt', 3)) +
    theme(axis.text.x = element_blank()) 

  PexcrDOC.p
  
  # N:P excretion
  NPexcrDOC.p <- plot_log_gam(gamNPDOC.pred) +
    geom_point(data = excr, aes(x = AmDOC, y = log10(massnorm.NP.excr),
                                   colour = Trophic.position), 
               size = point.size, alpha = point.alpha) +
    labs(x = expression(atop('DOC', paste('(mg C/L)'))),
         y = expression(atop(Log[10]~'mass-normalized', 
                             paste('N:P excretion (molar)')))) +
    scale_x_continuous(breaks = c(3, 5, 7, 9, 11))
  
  NPexcrDOC.p
  
  # ...DOM ----
  # N excretion
  NexcrDOM.p <- plot_gam(gamNDOM.pred) +
    geom_point(data = excr, aes(x = PC1, y = massnorm.N.excr,
                                   colour = Trophic.position), 
               size = point.size, alpha = point.alpha) +
    labs(x = '',
         y = '') +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank()) 
  NexcrDOM.p
  
  # P excretion
  PexcrDOM.p <- ggplot(excr, aes(x = PC1, y = massnorm.P.excr,
                                           colour = Trophic.position)) +
    geom_point(size = point.size, alpha = point.alpha) +
    geom_hline(data = excr.ss %>%
                 filter(Variable == 'massnorm.P.excr'), aes(yintercept = Mean),
                linetype = 'dashed', linewidth = line.size) +
    labs(x = '',
         y = '') +
    theme_classic(base_size = 10) +
    scale_color_manual(name = 'Trophic position',
                       labels = Trophic.labels,
                       values = met.brewer('Egypt', 3)) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank()) 

  PexcrDOM.p
  
  # N:P excretion
  NPexcrDOM.p <- plot_log_gam(gamNPDOM.pred) +
    geom_point(data = excr, aes(x = PC1, y = log10(massnorm.NP.excr),
                                   colour = Trophic.position), 
               size = point.size, alpha = point.alpha) +
    labs(x = expression(atop("DOM", paste("microbial-like"~~~~~~~~~~~~~~~~~~~
                                            "humic-like"))),
         y = '') +
    theme(axis.text.y = element_blank()) 
  
  NPexcrDOM.p
  
  # combine plots ----
  ggarrange(NexcrDOC.p, NexcrDOM.p, 
            PexcrDOC.p, PexcrDOM.p, 
            NPexcrDOC.p, NPexcrDOM.p,
            label.x = 0.25, label.y = 1,
            nrow = 3, ncol = 2, align = 'hv', heights = 1,
            common.legend = T,
            labels = c('(a)', '(b)', '(c)', '(d)', '(e)', '(f)'),
            font.label = list(size = 10))
  
  ggsave('figures/final-figures/Fig1.tiff', 
         width = 6, height = 8, 
         units = 'in', dpi = 600, compression = 'lzw')
  
  # Figure 2 ----
  # ...DOC excretion vs. DOC ----
  Cexcr.p <- plot_boxplot(excr.aov$massnorm.C.excr) +
    labs(x = '',
         y = expression(atop(Log[10]~'mass-normalized',
                             paste(DOC~excretion~(mg~C/g/h))))) +
    annotate("text", x = 2.5, y = 1, label = '*', size = text.size +.5) +
    geom_segment(x = 2, xend = 3, y = .9, yend = .9, linewidth = line.size) +
    geom_segment(x = 2, xend = 2, y = .9, yend = .8, linewidth = line.size) +
    geom_segment(x = 3, xend = 3, y = .9, yend = .8, linewidth = line.size) +
    theme(axis.text.x = element_blank())
  Cexcr.p
  
  # ...DOC:N excretion vs. DOC ----
  CNexcr.p <- plot_boxplot(excr.aov$massnorm.CN.excr) +
    labs(x = '',
         y = expression(atop(Log[10]~'DOC:N excretion',
                             paste((molar))))) +
    theme(axis.text.x = element_blank()) +
    annotate("text", x = 2.5, y = -1.15, label = '*', size = text.size + .5) +
    geom_segment(x = 2, xend = 3, y = -1.25, yend = -1.25, linewidth = line.size) +
    geom_segment(x = 2, xend = 2, y = -1.25, yend = -1.35, linewidth = line.size) +
    geom_segment(x = 3, xend = 3, y = -1.25, yend = -1.35, linewidth = line.size)
  CNexcr.p
  
  # ...DOC:P excretion vs. DOC ----
  CPexcr.p <- plot_boxplot(excr.aov$massnorm.CP.excr) +
    labs(x = 'DOC',
         y = expression(atop(Log[10]~'DOC:P excretion',
                             paste((molar))))) +
    theme(axis.title.x = element_text(vjust = -.4, hjust = 0.46))
  CPexcr.p
  
  # combine plots ----
  ggarrange(Cexcr.p, CNexcr.p, CPexcr.p,
            labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 10), label.x = 0.17, label.y = 1,
            nrow = 3, align = "v", legend = 'none')
  
  ggsave('figures/final-figures/Fig2.tiff', 
         width = 4, height = 6, 
         units = 'in', dpi = 600, compression = 'lzw')
  
  # Figure 3 ----
  DOMexcr.p <- ggplot(excr.DOM, aes(x = type, y = massnorm.excr)) +
    geom_boxplot(outlier.shape = NA, linewidth = line.size) +
    geom_jitter(aes(colour = DOC.level, fill = DOC.level), shape = 21, 
                size = point.size, stroke = .25, alpha = .6, width = .25) +
    labs(x = 'DOM optical parameters',
         y = 'Mass-normalized excretion rates') +
    geom_hline(aes(yintercept = 0), linetype = 'dashed') +
    guides(colour = guide_legend(override.aes = list(size = point.size + .5))) +
    scale_colour_manual(name = 'DOC',
                      labels = DOC.labels,
                      values = met.brewer('Greek', 3, direction = -1)) +
    scale_fill_manual(name = 'DOC',
                       labels = DOC.labels,
                       values = met.brewer('Greek', 3, direction = -1)) +
    theme_bw(base_size = 10) +
    theme(axis.title.x = element_text(vjust = -.4),
          axis.title.y = element_text(vjust = .7))
  DOMexcr.p
  
  ggsave('figures/final-figures/Fig3.tiff', 
         width = 7, height = 5, 
         units = 'in', dpi = 600, compression = 'lzw')
  
  # Figure 4 ----
  # ...low DOC ----
  nmds.l.p <- plot_nmds(nmds.l.scores) +
    stat_ellipse(level = .95, aes(colour = Trophic.position)) +
    scale_shape_manual(values = c(8, 16, 15, 17)) +
    scale_color_manual(values = c("#f5c34d", "steelblue4", "#43b284")) +
    scale_fill_manual(values = c("#f5c34d", "steelblue4", "#43b284")) +
    annotate("text", x = .65, y = 1.15,
             label = 'low DOC', colour = "#f5c34d", size = text.size) +
    annotate("text", x = -1.4, y = .45,
               label = 'omnivore', colour = "#43b284", size = text.size) +
    annotate("text", x = .4, y = .2,
               label = 'invertivore', colour = "steelblue4", size = text.size) +
    xlab('')
  nmds.l.p
  
  # ...medium DOC ----
  nmds.m.p <- plot_nmds(nmds.m.scores) +
    scale_shape_manual(values = c(8, 16)) +
    scale_color_manual(values = c("#8d1c06", "#dd5129")) +
    annotate("text", x = -1.6, y = .35, 
             label = 'medium DOC', colour = "#8d1c06", size = text.size) +
    annotate("text", x = .9, y = .1,
             label = 'invert/piscivore', colour = "#dd5129", size = text.size) +
    xlab('')
  nmds.m.p

  # ...high DOC ----
  nmds.h.p <- plot_nmds(nmds.h.scores) +
    scale_shape_manual(values = c(8, 16)) +
    scale_color_manual(values = c("#3c0d03", "#dd5129")) +
    annotate("text", x = 2, y = .55, 
             label = 'high DOC', colour = "#3c0d03", size = text.size) +
    annotate("text", x = -.2, y = -.6,
             label = 'invert/piscivore', colour = "#dd5129", size = text.size)
  nmds.h.p
  
  # combine plots ----
  ggarrange(nmds.l.p, nmds.m.p, nmds.h.p,
            labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 10), label.x = 0.13, label.y = 0.98,
            nrow = 3, align = "h", legend = 'none')
  ggsave('figures/final-figures/Fig4.tiff', 
         width = 10, height = 16, 
         units = 'cm', dpi = 600, compression = 'lzw')
  
  # Figure S1 ----
  # plot PCA
  PCA.p <- fviz_pca_biplot(pca,
                           repel = T,
                           label = 'var',
                           title = '',
                           ggtheme = theme_bw()) +
    xlim(-5,5) +
    ylim(-3,3) +
    labs(x = 'PC1 (75.7%)',
         y = 'PC2 (14.1%)') 
  PCA.p
  
  ggsave('figures/final-figures/FigS1.tiff', 
         width = 11, height = 11, 
         units = 'cm', dpi = 600, compression = 'lzw')
  