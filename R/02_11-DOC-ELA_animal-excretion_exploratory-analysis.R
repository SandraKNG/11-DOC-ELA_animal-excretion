  #### fish nutrient and DOC excretion across a lake DOC gradient ####
  # This code was created by S. Klemet-N'Guessan in 2022 and 2023
  # R version 4.3.0

  
  #### PRELIMINARY RESULTS ####
  # load libraries ----
  library(ggdist) # for stat_halfeye
  library(ggpubr) # for ggarrange function that allows to put all graphs on 1 page
  
  # labels
  Trophic.labels <- c('Invert/piscivore', 'Invertivore', 'Omnivore')
  Species.labels <- c('fathead minnow', 'pearl dace', 
                      'white sucker', 'yellow perch')
  
  # ..Trophic position ----
  # N excretion vs. DOC relative to trophic position
  Nexcrtp.p <- ggplot(excr.sp %>%  filter(Species.code != 'CTL1',
                                          Species.code != 'CTL2'), 
                      aes(x = AmBA, y = massnorm.N.excr.sp)) +
    geom_point(aes(color = Trophic.position), 
               size = 6) +   
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-normalized N excretion',
                             paste('(μg N/g/h)')))) +
    theme_grey(base_size = 26) +
    theme(axis.text = element_text(face = 'bold'),
          axis.line = element_line(size = 1),
          panel.grid = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'),
          legend.key.width = unit(3, 'lines'),
          legend.position = 'right') +
    scale_color_manual(name = 'Trophic position',
                       labels = Trophic.labels,
                       values = met.brewer('Egypt', 3))
  Nexcrtp.p
  
  # P excretion vs. DOC relative to trophic position
  Pexcrtp.p <- ggplot(excr.sp %>%  filter(Species.code != 'CTL1',
                                          Species.code != 'CTL2'), 
                      aes(x = AmDOC, y = massnorm.P.excr.sp)) +
    geom_point(aes(color = Trophic.position), 
               size = 6) +
    theme_grey(base_size = 26) +
    theme(axis.text = element_text(face = 'bold'),
          axis.line = element_line(size = 1),
          panel.grid = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'),
          legend.key.width = unit(3, 'lines'),
          legend.position = 'none') +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-normalized P excretion',
                             paste('(μg P/g/h)')))) +
    scale_color_manual(name = 'Trophic position',
                       labels = Trophic.labels,
                       values = met.brewer('Egypt', 3))
  Pexcrtp.p
  
  # N:P excretion
  NPexcrtp.p <- ggplot(excr.sp %>%  filter(Species.code != 'CTL1',
                                           Species.code != 'CTL2'), 
                       aes(x = AmDOC, y = massnorm.NP.excr.sp)) +
    geom_point(aes(color = Trophic.position), 
               size = 6) +
    theme_grey(base_size = 26) +
    theme(axis.text = element_text(face = 'bold'),
          axis.line = element_line(size = 1),
          panel.grid = element_blank(),
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(.15, .15, .15, .15, 'cm'),
          legend.key.height = unit(2, 'lines'),
          legend.key.width = unit(3, 'lines'),
          legend.position = 'none') +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-normalized',
                             paste('N:P excretion'~(μg~N/g/h))))) +
    scale_color_manual(name = 'Trophic position',
                       labels = Trophic.labels,
                       values = met.brewer('Egypt', 3))
  NPexcrtp.p
  
  
  # put them together
  ggarrange(Nexcrtp.p, Pexcrtp.p,
            labels = c("(a)", "(b)", "(c)"),
            font.label = list(size = 26), label.x = 0.26, label.y = 1,
            ncol = 2, align = "v", common.legend = T,
            legend = 'top')
  
  ggsave('figures/preliminary-figures/N_Pexcretion_tp.png', 
         width = 14, height = 7, 
         units = 'in', dpi = 600)
  
  # ..Species ----
  # N excretion vs. DOC relative to species
  Nexcrsp.p <- ggplot(excr.sp %>%  filter(Species.code != 'CTL1',
                                          Species.code != 'CTL2'), 
                      aes(x = AmDOC, y = massnorm.N.excr.sp)) +
    geom_point(aes(color = Species.code), 
               size = 4) +
    theme_classic(base_size = 20) +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-normalized',
                             paste(N~excretion~(μg~N/g/h))))) +
    scale_color_brewer(palette = "Set1",
                       name = 'Species',
                       labels = Species.labels)
  Nexcrsp.p
  
  # P excretion vs. DOC
  Pexcrsp.p <- ggplot(excr.sp %>%  filter(Species.code != 'CTL1',
                                          Species.code != 'CTL2'), 
                      aes(x = AmDOC, y = massnorm.P.excr.sp)) +
    geom_point(aes(color = Species.code), 
               size = 4) +
    theme_classic(base_size = 20) +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-normalized',
                             paste(P~excretion~(μg~P/g/h))))) +
    scale_color_brewer(palette = "Set1",
                       name = 'Species',
                       labels = Species.labels)
  Pexcrsp.p
  
  # put them together
  ggarrange(Nexcrsp.p, Pexcrsp.p,
            ncol = 2, align = "h", common.legend = T)
  
  ggsave('figures/preliminary-figures/N_Pexcretion_sp.png', 
         width = 14, height = 7, 
         units = 'in', dpi = 600)
  
  # ..fish average ----
  # N excretion vs.DOC
  Nexcrav.p <- ggplot(excr.smry, 
                      aes(x = AmDOC, y = massnorm.N.excr.sp)) +
    geom_point(size = 4) +
    theme_classic(base_size = 20) +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-normalized',
                             paste(N~excretion~(μg~N/g/h)))))
  Nexcrav.p
  
  # P excretion vs. DOC
  Pexcrav.p <- ggplot(excr.smry, 
                      aes(x = AmDOC, y = massnorm.P.excr.sp)) +
    geom_point(size = 4) +
    theme_classic(base_size = 20) +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-normalized',
                             paste(P~excretion~(μg~P/g/h))))) 
  Pexcrav.p
  
  # put them together
  ggarrange(Nexcrav.p, Pexcrav.p,
            ncol = 2, align = "h", common.legend = T)
  
  ggsave('preliminary results/N_Pexcretion_av.png', 
         width = 14, height = 7, 
         units = 'in', dpi = 600)

  
  # ..C excretion ----
  
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
  
  # ..DOM excretion vs ambient ----
  excr.nmds <- excr %>% select(c('Site.name', 'Species.code',
                                 'AmSUVA', 'AmSR', 'AmBA', 'AmFI', 'AmHIX', 'AmDOC',
                                 'SUVA.excretion', 'SR.excretion', 'BA.excretion',
                                 'FI.excretion', 'HIX.excretion'))
  
  excr.pcal <- excr %>% select(c('Site.name', 'Species.code',
                                 'AmSUVA', 'AmSR', 'AmBA', 'AmFI', 'AmHIX', 'AmDOC',
                                 'SUVA.excretion', 'SR.excretion', 'BA.excretion',
                                 'FI.excretion', 'HIX.excretion'))
  
  excr.pcal <- excr.pcal %>% rename(SUVA254 = AmSUVA,
                                    βα = AmBA,
                                    SR = AmSR,
                                    FI = AmFI,
                                    HIX = AmHIX,
                                    DOC = AmDOC)
  
  excr.pca <- excr.pcal %>% group_by(Site.name) %>%
    summarise(SUVA254 = mean(SUVA254),
              βα = mean(βα),
              SR = mean(SR),
              FI = mean(FI),
              HIX = mean(HIX),
              DOC = mean(DOC))
  # excr.nmds <- excr.nmds %>% filter(!is.na(SUVA.excretion)) %>% 
  #   rename(SUVA254 = AmSUVA,
  #          βα = AmBA,
  #          SR = AmSR,
  #          FI = AmFI,
  #          HIX = AmHIX,
  #          DOC = AmDOC)
  
  excr.nmdsm <- as.matrix(excr.nmds)
  
  pca <- princomp(excr.pca[, 2:7], cor = TRUE, scores = TRUE)
  biplot(pca)
  # Extract PC1 scores
  pca.PC1 <- pca$scores[,1]
  
  # plot PCA1
  library(ggbiplot) # to plot PCA results, only load it when need it for plotting
  PCA.p <- ggbiplot(pca, obs.scale = 1, var.scale = 1) +
    labs(x = 'PC1 (70.4%)',
         y = 'PC2 (19.6%)') +
    # xlim(-4,4.5) + ylim(-4.5,5) +
    theme_grey(base_size = 13) +
    theme(panel.grid = element_blank())
  
  PCA.p
  
  ggsave('figures/preliminary-figures/PCA_DOM.png', 
         width = 11, height = 9, 
         units = 'cm', dpi = 600)
  
  detach("package:ggbiplot")
  
  nmds <- metaMDS(excr.nmds[, 3:13], distance = 'bray')
  plot(nmds)
  
  # extract NMDS scores (x and y coordinates)
  data.scores = as.data.frame(scores(nmds))
  data.scores = cbind(excr.nmds %>% filter(!is.na(SUVA.excretion)) %>% 
                        select(Site.name, Species.code))
  data.scores <- mutate(Site.name = excr.nmds$Site.name)
