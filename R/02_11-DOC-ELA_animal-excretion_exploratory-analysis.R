  #### fish nutrient and DOC excretion across a lake DOC gradient ####
  # This code was created by S. Klemet-N'Guessan in 2022 and 2023
  # R version 4.3.0

  
  #### PRELIMINARY RESULTS ####
  # load libraries ----
  library(ggdist) # for stat_halfeye
  library(ggpubr) # for ggarrange function that allows to put all graphs on 1 page
  library(MetBrewer)
  
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
  
  # DOC vs DOC levels
  ggboxplot(excr %>%  dplyr::filter(!is.na(massnorm.C.excr)),
            "DOC.level", "massnorm.C.excr") 
  excr %>% 
    select(c(C.excretion.rate, massnorm.C.excr, DOC.level)) %>% 
    group_by(DOC.level) %>% 
    identify_outliers(C.excretion.rate)
  
  # normality check
  # Build the linear model
  model  <- lm(log10(massnorm.C.excr) ~ DOC.level, data = excr.aov)
  # Create a QQ plot of residuals
  ggqqplot(residuals(model))
  
  # Compute Shapiro-Wilk test of normality
  shapiro_test(residuals(model))
  
  # QQ plot for each group level
  ggqqplot(excr.aov, 'massnorm.C.excr', facet.by = 'DOC.level')
  
  # check homogeneity of variances
  excr.aov %>% levene_test(massnorm.C.excr ~ DOC.level)


  
  
  
  