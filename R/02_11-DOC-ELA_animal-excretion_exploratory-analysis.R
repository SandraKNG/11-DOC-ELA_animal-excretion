  #### fish nutrient and DOC excretion across a lake DOC gradient ####
  # This code was created by S. Klemet-N'Guessan in 2022 and 2023
  # R version 4.3.0

  
  #### PRELIMINARY RESULTS ####
  # load libraries ----
  library(ggdist) # for stat_halfeye
  library(ggpubr) # for ggarrange function that allows to put all graphs on 1 page
  library(vegan) # for NMDS
  
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
  
  # DOM excretion only ----
  ggplot(excr.DOM, aes(x = type, y = massnorm.excr)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha = .3, width = .25) +
    theme_light()
  
  excr.var <- excr.var %>% mutate(across(where(is.numeric), 
                                         ~ if_else(. < 0, 0, .)))
  t.test(excr.var$massnorm.SR.excr, mu = 0, alternative = c("greater"))
  
  # ..DOM excretion vs ambient ----
  # function to rename some DOM characteristics
  rename_DOM <- function(df) {
    df <- df %>% 
      rename(
      SUVA254 = AmSUVA,
      βα = AmBA,
      SR = AmSR,
      FI = AmFI,
      HIX = AmHIX,
      DOC = AmDOC,
      C1 = AmC1,
      C2 = AmC2,
      C3 = AmC3,
      C4 = AmC4,
      C5 = AmC5,
      C7 = AmC7,
    )
    return(df)
  }
  
  # PCA ----
  # prepare PCA dataset
  excr.pca <- excr %>% group_by(Site.name) %>%
    summarise(across(c(
      starts_with('Am'),-AmC6,-starts_with(c(
        'AmA', 'AmP', 'AmS2', 'Ambi', 'AmS3', 'AmHis', 'AmR'
      ))
    ),
    \(x) mean(x, na.rm = TRUE))) %>% 
    rename_DOM()
  # do PCA
  # perhaps should use %?
  pca <- princomp(excr.pca[, 3:13], cor = TRUE, scores = TRUE)
  biplot(pca)
  # Extract PC1 scores
  pca.PC1 <- pca$scores[,1]
  
  # plot PCA1
  library(ggbiplot) # to plot PCA results, only load it when need it for plotting
  PCA.p <- ggbiplot(pca, obs.scale = 1, var.scale = 1) +
    # labs(x = 'PC1 (70.4%)',
    #      y = 'PC2 (19.6%)') +
    # xlim(-4,4.5) + ylim(-4.5,5) +
    theme_bw(base_size = 13) +
    theme(panel.grid = element_blank())
  
  PCA.p
  
  ggsave('figures/preliminary-figures/PCA_DOM.png', 
         width = 11, height = 9, 
         units = 'cm', dpi = 600)
  
  detach("package:ggbiplot")
  
  # NMDS ----
  # prepare NMDS dataset
  excr.amb <- excr.pca %>% 
    rename(ID = Site.name) %>% 
    filter(ID %in% c('L222', 'L224', 'L239')) %>% 
    dplyr::mutate(Site.name = c('L222', 'L224', 'L239'),
                  Species.code = c('L222', 'L224', 'L239'))
  excr.nmds <- excr %>% 
    #group_by(Site.name, Species.code) %>%
    select(c(
      ID, Site.name, Species.code,
      ends_with('excretion.rate'),
      -starts_with(c('N.e', 'P.e'))
      # starts_with('Am'),
      # -AmC6,
      # -starts_with(
      #   c('AmA', 'AmP', 'AmS2', 'Ambi', 'N.e',
      #     'P.e', 'AmS3', 'AmHis', 'AmR')
      )
    ) %>% 
  dplyr::filter(
    !between(ID, 935, 946),
    !ID %in% c(1008, 1012, 1013, 115, 901, 904, 907, 911, 
               907, 927, 928, 929, 930)
  ) %>% 
    mutate(across(where(is.numeric),
                  ~ if_else(. < 0, 0, .)),
           ID = as.character(ID)) %>% 
    rename_with(~ sub(".excretion.rate", "", .), .cols = where(is.numeric)) %>% 
    rename(
      DOC = C,
      βα = BA,
      SUVA254 = SUVA
    ) %>% 
    dplyr::filter(
      !is.na(DOC)) %>% 
    bind_rows(excr.amb)
    # summarise(across(everything(), mean, na.rm = TRUE)) %>%
    # rename_DOM() %>%
    # dplyr::filter(!is.na(C.excretion.rate)) %>% 
    # pivot_longer(
    #   cols = C.excretion.rate:C7.excretion.rate,
    #   #names_pattern = '(*.)excretion.rate',
    #   names_to = 'measurement',
    #   values_to = 'value'
    # ) %>%
    # pivot_wider(
    #   id_cols = c('Site.name', 'Species.code'),
    #   names_from = ID,
    #   values_from = value)
  
  # transform NMDS dataset to a matrix
  excr.nmds.m <- excr.nmds %>% 
    select(-c(ID, Species.code, Site.name)) %>% 
    as.matrix() 
  
  # do NMDS
  set.seed(1)
  nmds <- metaMDS(excr.nmds.m, distance = 'bray', k = 2, trymax = 100)
  stressplot(nmds)
  
  # First create a data frame of the scores from the individual sites.
  # This data frame will contain x and y values for where sites are located.
  nmds.scores <- scores(nmds)$sites %>% 
    as.tibble(rownames = 'sample') %>% 
    dplyr::mutate(
      ID = excr.nmds$ID,
      Site.name = excr.nmds$Site.name,
      Species.code = excr.nmds$Species.code
    )
  
  # NMDS plot
  ggplot(nmds.scores, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(colour = Species.code), size = 4) +
    scale_color_brewer(palette = 'Set1') +
    theme_bw()

  