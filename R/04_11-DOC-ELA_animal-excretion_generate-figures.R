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
  
  pred_df <- function(olddf, x) {
    if (any(x == 'excr$AmDOC')) {
      newdf <- with(olddf, 
                    expand.grid(x = c(seq(min(x), max(x),
                                          length = 100), 
                                      rep(median(x), 100)),
                                Trophic.position = levels(Trophic.position)))
    } else {
      newdf <- with(olddf, 
                    expand.grid(x = c(seq(min(x), max(x),
                                          length = 100), 
                                      rep(median(x), 100)),
                                Trophic.position = levels(Trophic.position)))
    }
    function(mod) {newdf <- cbind(newdf, 
                                  predict(mod,newdf,
                                          se.fit = TRUE,
                                          type = "link"))
    ilink <- inv_link(mod)
    newdf <- newdf %>% mutate(response = ilink(fit),
                              upper = ilink(fit + (1.96 * se.fit)),
                              lower = ilink(fit - (1.96 * se.fit)))
    return(newdf)
    }
  }
  
  # Prediction function for vertebrate classification with  pred on log scale 
  pred_log_df <- function(olddf, x) {
    if (any(x == 'excr$AmDOC')) {
      newdf <- with(olddf, 
                    expand.grid(x = c(seq(min(x), max(x),
                                          length = 100), 
                                      rep(median(x), 100)),
                                group = levels(group)))
    } else {
      newdf <- with(olddf, 
                    expand.grid(x = c(seq(min(x), max(x),
                                          length = 100), 
                                      rep(median(x), 100)),
                                group = levels(group)))
    }
    
    function(mod) {newdf <- cbind(newdf, 
                                  predict(mod,newdf,
                                          se.fit = TRUE,
                                          type = "response"))
    newdf <- newdf %>% mutate(upper = fit + (1.96 * se.fit),
                              lower = fit - (1.96 * se.fit))
    return(newdf)
    }
  }
  
  # Storing prediction functions with different inputs
  pred.DOC <- pred_df(excr, excr$AmDOC)
  pred.DOM <- pred_df(excr, excr$PC1)
  
  # *Prediction models ----
  gamNDOC.pred <- pred.DOC(gamNDOC)
  
  
  # set up plotting parameters ----
  Trophic.labels <- c('Invert/piscivore', 'Invertivore', 'Omnivore')
  Species.labels <- c('fathead minnow', 'pearl dace', 
                      'white sucker', 'yellow perch')
  point.size = 1.5
  line.size = 0.5
  
  # Figure 1 ----
  # ..Option 1 ----
  # N excretion
  NexcrDOC.p <- ggplot(gamNDOC.pred, 
                         aes(x = AmDOC, y = fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = .2) +
    geom_line(linewidth = 1.5) +
    geom_point(data = excr.sp, aes(x = AmDOC, y = massnorm.N.excr.sp,
                                   colour = Trophic.position), 
               size = 5) +
    labs(x = '',
         y = expression(atop('Mass-normalized',
                             paste(N~excretion~(μg~N/g/h))))) +
    theme_grey(base_size = 22) +
    theme(axis.text = element_text(face = 'bold'),
          text = element_text(family = "Helvetica"),
          axis.line = element_line(size = 0.5),
          panel.grid = element_blank(),
          plot.margin = unit(c(-0.1,0.7,0,0), "lines"),
          axis.text.x = element_blank(),
          strip.text.y = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .20, .15, 'cm'),
          legend.key.height = unit(1, 'lines'), 
          legend.key.width = unit(2, 'lines'),
          legend.position = 'top') +
    scale_color_manual(name = 'Trophic position',
                       labels = Trophic.labels,
                       values = met.brewer('Egypt', 3)) 
  
  
  NexcrDOC.p
  
  # P excretion
  PexcrDOC.p <- ggplot(gamPDOC.pred, 
                       aes(x = AmDOC, y = fit)) +
    # geom_ribbon(aes(ymin = lower, ymax = upper),
    #             alpha = .1) +
    # geom_line(linewidth = 0.5) +
    geom_point(data = excr.sp, aes(x = AmDOC, y = massnorm.P.excr.sp,
                                   colour = Trophic.position), 
               size = 5) +
    geom_hline(data = excr.ss %>%
                 filter(Variable == 'massnorm.P.excr'), 
               aes(yintercept = Mean), linetype = 'dashed', size = 1.5) +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-normalized',
                             paste(P~excretion~(μg~P/g/h))))) +
    theme_grey(base_size = 22) +
    theme(axis.text = element_text(face = 'bold'),
          text = element_text(family = "Helvetica"),
          axis.line = element_line(size = 0.5),
          panel.grid = element_blank(),
          plot.margin = unit(c(-0.1,0.7,0,0), "lines"),
          #axis.text.x = element_blank(),
          strip.text.y = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .20, .15, 'cm'),
          legend.key.height = unit(1, 'lines'), 
          legend.key.width = unit(2, 'lines'),
          legend.position = 'top') +
    scale_color_manual(name = 'Trophic position',
                       labels = Trophic.labels,
                       values = met.brewer('Egypt', 3)) 
  
  PexcrDOC.p
  
  # P excretion
  NPexcrDOC.p <- ggplot(gamNPDOC.pred, 
                       aes(x = AmDOC, y = fit)) +
    # geom_ribbon(aes(ymin = lower, ymax = upper),
    #             alpha = .1) +
    # geom_line(linewidth = 0.5) +
    geom_point(data = excr.sp, aes(x = AmDOC, y = massnorm.NP.excr.sp), 
               size = 5, alpha = .5) +
    geom_hline(data = excr.ss %>%
                 filter(Variable == 'massnorm.NP.excr'), 
               aes(yintercept = Mean), linetype = 'dashed', size = 1.5) +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop(Log[10]~'mass-normalized',
                             paste(N:P~excretion~(molar))))) +
    theme_grey(base_size = 24) +
    theme(axis.text = element_text(face = 'bold'),
          text = element_text(family = "Helvetica"),
          axis.line = element_line(size = 0.5),
          panel.grid = element_blank(),
          plot.margin = unit(c(-0.1,0.7,0,0), "lines"),
          # axis.text.x = element_blank(),
          strip.text.y = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .20, .15, 'cm'),
          legend.key.height = unit(1, 'lines'), 
          legend.key.width = unit(2, 'lines'),
          legend.position = 'top') 
  
  NPexcrDOC.p
  
  # put them together
  ggarrange(NexcrDOC.p, PexcrDOC.p, 
            labels = c("(a)", "(b)"),
            font.label = list(size = 22), label.x = 0.22, label.y = 1,
            nrow = 2, align = "hv", common.legend = T,
            legend = 'right')
  
  ggsave('figures/preliminary-figures/N_Pexcretion_gam_tp.png', 
         width = 10, height = 10, 
         units = 'in', dpi = 300)
  
  ggarrange(NexcrDOC.p, PexcrDOC.p, NPexcrDOC.p,
            labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 26), label.x = 0.22, label.y = 1,
            nrow = 3, align = "hv", common.legend = T)
  
  ggsave('figures/preliminary-figures/N_Pexcretion_gam2.png', 
         width = 7, height = 14, 
         units = 'in', dpi = 300)
  
  # ..Option 2 ----
  # N excretion
  NexcrDOCtp.p <- ggplot(gamNDOCtp.pred, 
                          aes(x = AmDOC, y = massnorm.N.excr, 
                              group = Trophic.position)) +
    # geom_ribbon(aes(ymin = lower, ymax = upper, fill = Trophic.position),
    #             alpha = .1) +
    geom_line(aes(colour =Trophic.position), linewidth = 0.5) +
    geom_point(data = excr, aes(x = AmDOC, y = massnorm.N.excr,
                                color = Trophic.position), 
               size = 5, alpha = .4) +
    # # geom_hline(data = excr.ss %>%  
    #              filter(Variable == 'Log10.massnorm.N.excr',
    #                     .group == 'Vert.class=Vert'), aes(yintercept = Mean), 
    #            color = "#377EB8", linetype = 'dashed', size = 0.5) +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-normalized',
                             paste(N~excretion~(μg~N/g/h))))) +
    theme_classic(base_size = 26) +
    theme(axis.text = element_text(face = 'bold'),
      text = element_text(family = "Helvetica"),
      axis.line = element_line(size = 0.5),
      panel.grid = element_blank(),
      plot.margin = unit(c(-0.1,0.7,0,0), "lines"),
      axis.text.x = element_blank(),
      strip.text.y = element_blank(),
      legend.title = element_text(face = 'bold'),
      legend.margin = margin(.15, .15, .20, .15, 'cm'),
      legend.key.height = unit(1, 'lines'), 
      legend.key.width = unit(2, 'lines'),
      legend.position = 'top') +
    scale_color_manual(name = 'Trophic position',
                       labels = Trophic.labels,
                       values = met.brewer('Egypt', 3)) +
    scale_fill_manual(name = 'Trophic position',
                       labels = Trophic.labels,
                       values = met.brewer('Egypt', 3))
  
  NexcrDOCtp.p
  
  # N excretion
  PexcrDOCtp.p <- ggplot(gamPDOCtp.pred, 
                         aes(x = AmDOC, y = massnorm.P.excr, 
                             group = Trophic.position)) +
    # geom_ribbon(aes(ymin = lower, ymax = upper, fill = Trophic.position),
    #             alpha = .1) +
    geom_line(aes(colour =Trophic.position), linewidth = 0.5) +
    geom_point(data = excr, aes(x = AmDOC, y = massnorm.P.excr,
                                color = Trophic.position), 
               size = 5, alpha = .4) +
    # # geom_hline(data = excr.ss %>%  
    #              filter(Variable == 'Log10.massnorm.N.excr',
    #                     .group == 'Vert.class=Vert'), aes(yintercept = Mean), 
    #            color = "#377EB8", linetype = 'dashed', size = 0.5) +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-normalized',
                             paste(P~excretion~(μg~P/g/h))))) +
    theme_classic(base_size = 26) +
    theme(axis.text = element_text(face = 'bold'),
          text = element_text(family = "Helvetica"),
          axis.line = element_line(size = 0.5),
          panel.grid = element_blank(),
          plot.margin = unit(c(-0.1,0.7,0,0), "lines"),
          axis.text.x = element_blank(),
          strip.text.y = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .20, .15, 'cm'),
          legend.key.height = unit(1, 'lines'), 
          legend.key.width = unit(2, 'lines'),
          legend.position = 'top') +
    scale_color_manual(name = 'Trophic position',
                       labels = Trophic.labels,
                       values = met.brewer('Egypt', 3)) +
    scale_fill_manual(name = 'Trophic position',
                      labels = Trophic.labels,
                      values = met.brewer('Egypt', 3))
  
  PexcrDOCtp.p
  
  # Figure 2 ----
  # C excretion vs. DOC
  Cexcr.p <- ggplot(excr %>%  filter(!is.na(DOC.level),
                                     C.excretion.rate > 0), 
                    aes(x = DOC.level, y = log10(massnorm.C.excr), 
                        group = DOC.level)) +
    stat_halfeye(adjust = .5, width = .7, .width = 0,justification = -.4,
                 alpha = .3, aes(fill = DOC.level)) + 
    geom_boxplot(width = .4, size = 1.1, outlier.shape = NA, 
                 aes(fill = DOC.level, colour = DOC.level), alpha = .2) +
    geom_point(size = 5, alpha = .3, 
               position = position_jitter(seed = 1, width = .1),
               aes(colour = DOC.level)) +
    labs(x = '',
         y = expression(atop(Log[10]~'mass-normalized',
                             paste(C~excretion~(mg~C/g/h))))) +
    scale_x_discrete(labels = c('low', 'medium', 'high')) +
    theme_grey(base_size = 26) +
    # theme(axis.text = element_text(face = 'bold'),
    #       axis.line = element_line(size = 1),
    #       panel.grid = element_blank(),
    #       axis.text.x = element_blank(),
    #       legend.title = element_text(face = 'bold'),
    #       legend.margin = margin(.15, .15, .15, .15, 'cm'),
    #       legend.key.height = unit(2, 'lines'),
    #       legend.key.width = unit(3, 'lines'),
    #       legend.position = 'none') +
    scale_fill_manual(name = 'DOC level',
                      labels = c('low', 'medium', 'high'),
                      values = met.brewer('Greek', 3, direction = -1)) +
    scale_colour_manual(name = 'DOC level',
                        labels = c('low', 'medium', 'high'),
                        values = met.brewer('Greek', 3, direction = -1)) 
  Cexcr.p
  
  # C:N excretion vs. DOC
  CNexcr.p <- ggplot(excr %>%  filter(!is.na(DOC.level),
                                      C.excretion.rate > 0), 
                     aes(x = DOC.level, y = log10(massnorm.CN.excr), 
                         group = DOC.level)) +
    stat_halfeye(adjust = .5, width = .7, .width = 0,justification = -.4,
                 alpha = .3, aes(fill = DOC.level)) + 
    geom_boxplot(width = .4, size = 1.1, outlier.shape = NA, 
                 aes(fill = DOC.level, colour = DOC.level), alpha = .2) +
    geom_point(size = 5, alpha = .3, 
               position = position_jitter(seed = 1, width = .1),
               aes(colour = DOC.level)) +
    labs(x = '',
         y = expression(atop(Log[10]~'C:N excretion',
                             paste((molar))))) +
    scale_x_discrete(labels = c('low', 'medium', 'high')) +
    theme_grey(base_size = 26) +
    theme(axis.text = element_text(face = 'bold'),
          axis.line = element_line(size = 1),
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'),
          legend.key.width = unit(3, 'lines'),
          legend.position = 'none') +
    scale_fill_manual(name = 'DOC level',
                      labels = c('low', 'medium', 'high'),
                      values = met.brewer('Greek', 3, direction = -1)) +
    scale_colour_manual(name = 'DOC level',
                        labels = c('low', 'medium', 'high'),
                        values = met.brewer('Greek', 3, direction = -1)) 
  CNexcr.p
  
  # C:P excretion vs. DOC
  CPexcr.p <- ggplot(excr %>%  filter(!is.na(DOC.level),
                                      C.excretion.rate > 0), 
                     aes(x = DOC.level, y = log10(massnorm.CP.excr), 
                         group = DOC.level)) +
    stat_halfeye(adjust = .5, width = .7, .width = 0,justification = -.4,
                 alpha = .3, aes(fill = DOC.level)) + 
    geom_boxplot(width = .4, size = 1.1, outlier.shape = NA, 
                 aes(fill = DOC.level, colour = DOC.level), alpha = .2) +
    geom_point(size = 4, alpha = .3, 
               position = position_jitter(seed = 1, width = .1),
               aes(colour = DOC.level)) +
    labs(x = 'DOC',
         y = expression(atop(Log[10]~'C:P excretion',
                             paste((molar))))) +
    scale_x_discrete(labels = c('low', 'medium', 'high')) +
    theme_grey(base_size = 26) +
    theme(axis.text = element_text(face = 'bold'),
          axis.line = element_line(size = 1),
          panel.grid = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'),
          legend.key.width = unit(3, 'lines'),
          legend.position = 'none') +
    scale_fill_manual(name = 'DOC level',
                      labels = c('low', 'medium', 'high'),
                      values = met.brewer('Greek', 3, direction = -1)) +
    scale_colour_manual(name = 'DOC level',
                        labels = c('low', 'medium', 'high'),
                        values = met.brewer('Greek', 3, direction = -1)) 
  CPexcr.p
  
  
  # put them together
  ggarrange(Cexcr.p, CNexcr.p, CPexcr.p,
            labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 26), label.x = 0.22, label.y = 1,
            nrow = 3, align = "v", legend = 'none')
  
  ggsave('figures/preliminary-figures/C_N_Pexcretion.png', 
         width = 7, height = 14, 
         units = 'in', dpi = 600)
  
  # Figure 3 ----
  DOMexcr.p <- ggplot(excr.DOM, aes(x = type, y = massnorm.excr)) +
    geom_boxplot(outlier.shape = NA, linewidth = line.size) +
    geom_jitter(aes(colour = DOC.level, fill = DOC.level), shape = 21, 
                size = point.size, stroke = .25, alpha = .6, width = .25) +
    labs(x = '',
         y = 'Mass-normalized excretion rates') +
    geom_hline(aes(yintercept = 0), linetype = 'dashed') +
    theme(legend.key.height= unit(1, 'lines'),
          legend.key.width= unit(2, 'lines')) +
    # scale_y_continuous(trans='log10',
    #                    breaks=trans_breaks('log10', function(x) 10^x),
    #                    labels=trans_format('log10', math_format(10^.x))) +
    scale_colour_manual(name = 'DOC',
                      labels = c('low', 'medium', 'high'),
                      values = met.brewer('Greek', 3, direction = -1)) +
    scale_fill_manual(name = 'DOC',
                       labels = c('low', 'medium', 'high'),
                       values = met.brewer('Greek', 3, direction = -1)) +
    theme_bw() 
  DOMexcr.p
  
  ggsave('figures/final-figures/DOM_charact_excretion.png', 
         width = 7, height = 5, 
         units = 'in', dpi = 600)
  
  t.test(excr.var$massnorm.SR.excr, mu = 0, alternative = c("greater"))
  
  # Figure 4 ----
  # NMDS plot
  ggplot(nmds.scores, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(colour = Species.code), size = 2) +
    # scale_color_brewer(palette = 'Set1') +
    scale_color_manual(values = met.brewer("Archambault")) +
    theme_bw()
  
  ggsave('figures/preliminary-figures/NMDS.png', 
         width = 11, height = 7, 
         units = 'cm', dpi = 600)
  
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
  
  ggsave('figures/final-figures/PCA_DOM.png', 
         width = 11, height = 11, 
         units = 'cm', dpi = 600)
  