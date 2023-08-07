  #### fish nutrient and DOC excretion across a lake DOC gradient ####
  # This code was created by S. Klemet-N'Guessan in 2022 and 2023
  # R version 4.3.0
  
  ################################ FINAL FIGURES #################################
  # load libraries ----
  library(ggdist) # for stat_halfeye
  library(ggpubr) # for ggarrange function that allows to put all graphs on 1 page
  
  # Set up prediction data ----
  # without trophic position and with population averages
  gamNDOC.pred <- with(excr, 
                        expand.grid(AmDOC = c(seq(min(AmDOC), max(AmDOC),
                                                  length = 100), 
                                              rep(median(AmDOC), 100))))
  
  gamPDOC.pred <- with(excr, 
                       expand.grid(AmDOC = c(seq(min(AmDOC), max(AmDOC),
                                                 length = 100), 
                                             rep(median(AmDOC), 100))))
  
  gamNPDOC.pred <- with(excr, 
                       expand.grid(AmDOC = c(seq(min(AmDOC), max(AmDOC),
                                                 length = 100), 
                                             rep(median(AmDOC), 100))))
  
  # with trophic position
  gamNDOCtp.pred <- with(excr, 
                       expand.grid(AmDOC = c(seq(min(AmDOC), max(AmDOC),
                                                 length = 100), 
                                             rep(median(AmDOC), 100)),
                                   Trophic.position = levels(Trophic.position)))
  
  gamPDOCtp.pred <- with(excr, 
                       expand.grid(AmDOC = c(seq(min(AmDOC), max(AmDOC),
                                                 length = 100), 
                                             rep(median(AmDOC), 100)),
                                   Trophic.position = levels(Trophic.position)))
  
  # Model with DOC and DOM separately ----
  # without trophic position and with population averages
  gamNDOC.pred <- cbind(gamNDOC.pred,
                         predict(gamNDOC, 
                                 gamNDOC.pred, 
                                 se.fit = TRUE, 
                                 type = "response"))
  
  gamNDOC.pred <- gamNDOC.pred %>% mutate(upper = fit + (1.96 * se.fit),
                                            lower = fit - (1.96 * se.fit))
  
  gamPDOC.pred <- cbind(gamPDOC.pred,
                        predict(gamPDOC, 
                                gamPDOC.pred, 
                                se.fit = TRUE, 
                                type = "response"))
  
  gamPDOC.pred <- gamPDOC.pred %>% mutate(upper = fit + (1.96 * se.fit),
                                            lower = fit - (1.96 * se.fit))
  
  gamNPDOC.pred <- cbind(gamNPDOC.pred,
                        predict(gamNPDOC, 
                                gamNPDOC.pred, 
                                se.fit = TRUE, 
                                type = "response"))
  
  gamNPDOC.pred <- gamNPDOC.pred %>% mutate(upper = fit + (1.96 * se.fit),
                                          lower = fit - (1.96 * se.fit))
  
  # with trophic position
  gamNDOCtp.pred <- cbind(gamNDOCtp.pred,
                        predict(gamNDOC.tp, 
                                gamNDOCtp.pred, 
                                se.fit = TRUE, 
                                type = "link"))
  ilink <- inv_link(gamNDOC.tp)
  gamNDOCtp.pred <- gamNDOCtp.pred %>% mutate(massnorm.N.excr = ilink(fit),
                                          upper = ilink(fit + (1.96 * se.fit)),
                                          lower = ilink(fit - (1.96 * se.fit)))
  
  gamPDOCtp.pred <- cbind(gamPDOCtp.pred,
                        predict(gamPDOC.tp, 
                                gamPDOCtp.pred, 
                                se.fit = TRUE, 
                                type = "link"))
  ilink <- inv_link(gamPDOC.tp)
  gamPDOCtp.pred <- gamPDOCtp.pred %>% mutate(massnorm.P.excr = ilink(fit),
                                          upper = ilink(fit + (1.96 * se.fit)),
                                          lower = ilink(fit - (1.96 * se.fit)))
  
  # Figure 1 ----
  # labels
  Trophic.labels <- c('Invert/piscivore', 'Invertivore', 'Omnivore')
  Species.labels <- c('fathead minnow', 'pearl dace', 
                      'white sucker', 'yellow perch')
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
  