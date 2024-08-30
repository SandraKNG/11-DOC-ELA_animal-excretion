  #### Fish supply distinct nutrients and dissolved organic matter composition ####
  ### relative to ambient concentrations in northern lakes ####
  
  # This code was created by S. Klemet-N'Guessan in 2022-2024
  # R version 4.3.0

  
  #### PRELIMINARY RESULTS ####
  # load libraries ----
  library(ggdist) # for stat_halfeye
  library(ggpubr) # for ggarrange function that allows to put all graphs on 1 page
  library(MetBrewer)
  library(performance)
  
  # labels
  Trophic.labels <- c('Invert/piscivore', 'Invertivore', 'Omnivore')
  Species.labels <- c('fathead minnow', 'pearl dace', 
                      'white sucker', 'yellow perch')
  
  
  # Look at N/P excretion vs. mass ----
  #..Fathead minnows only  ----
  because coeffs are too high for N excr (>1) + too low for P excr (<0.3)
  N excretion
  ggplot(excr %>%  filter(Species.code == 'FM',
                          Log10.N.excretion.rate > 0.5),
         aes(x = Log10.mass, y = Log10.N.excretion.rate)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    theme(legend.position = 'none')
  modelN.FM <- lm(Log10.N.excretion.rate ~ Log10.mass,
               data = excr %>%
                 filter(Species.code == 'FM',
                        Log10.N.excretion.rate > 0.5))
  modelN.FM$coefficients["Log10.mass"]

  # P excretion
  ggplot(excr %>%  filter(Species.code == 'FM'),
         aes(x = Log10.mass, y = Log10.P.excretion.rate)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    theme(legend.position = 'none')
  modelP.FM <- lm(Log10.P.excretion.rate ~ Log10.mass,
                  data = excr %>%
                    filter(Species.code == 'FM'))
  modelP.FM$coefficients["Log10.mass"]
  
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
  
  # testing things
  
  ggplot(excr.vol, aes(x = vol.Cexcr, y = surf.Cexcr_d)) +
    geom_point(size = 2) +
    theme_classic() +
    labs(x = 'Volumetric C excretion (mg/L)',
         y = 'Surface C excretion rate (mg/d)') 
  
  ggplot(excr.vol, aes(x = AmDOC, y = turnover.N.time_yr)) +
    geom_point(size = 2) +
    geom_segment( aes(x=AmDOC, xend=AmDOC, y=0, yend=turnover.N.time_yr)) +
    # geom_bar(stat = "identity") +
    theme_bw() +
    scale_y_log10() +
    labs(x = 'DOC (mg/L)',
         y = 'N turnover time (year)')
  
  ggplot(excr.vol, aes(x = AmDOC, y = turnover.P.time_d)) +
    geom_point(size = 2) +
    geom_segment( aes(x=AmDOC, xend=AmDOC, y=0, yend=turnover.P.time_d)) +
    # geom_bar(stat = "identity") +
    theme_bw() +
    scale_y_log10() +
    labs(x = 'DOC (mg/L)',
         y = 'P turnover time (days)')
  
  ggplot(excr.vol, aes(x = AmDOC, y = turnover.C.time_yr)) +
    geom_point(size = 2) +
    geom_segment( aes(x=AmDOC, xend=AmDOC, y=0, yend=turnover.C.time_yr)) +
    # geom_bar(stat = "identity") +
    theme_bw() +
    scale_y_log10() +
    labs(x = 'DOC (mg/L)',
         y = 'C turnover time (year)')
  
  # Do excretion rates change with incubation temperature (temperature effect)?
  
  ggplot(excr, aes(x = Temp, y = log10(N.excretion.rate))) +
    geom_point(size = 2) +
    geom_smooth() +
    theme_classic() 
  
  ggplot(excr, aes(x = Temp, y = P.excretion.rate)) +
    geom_point(size = 2) +
    geom_smooth() +
    theme_classic() 
  
  ggplot(excr, aes(x = Temp, y = log10(massnorm.N.excr))) +
    geom_point(size = 2) +
    geom_smooth(method = lm) +
    theme_classic()
  
  ggplot(excr, aes(x = Temp, y = log10(massnorm.P.excr))) +
    geom_point(size = 2) +
    geom_smooth(method = lm) +
    theme_classic()
  
  temp.N <- lm(log10(massnorm.N.excr) ~ Temp, data = excr)
  anova(temp.N)
  temp.N$coefficients
  check_model(temp.N)
  
  temp.P <- lm(log10(massnorm.P.excr) ~ Temp, data = excr)
  anova(temp.P)
  temp.P$coefficients
  check_model(temp.P)
  
  #  Do excretion rates change with incubation time (stress effect)?
  ggplot(excr, aes(x = Time.elapsed, y = log10(massnorm.N.excr))) +
    geom_point(size = 2) +
    geom_smooth(method = lm) +
    theme_classic()
  
  time.N <- lm(log10(massnorm.N.excr) ~ Time.elapsed, data = excr)
  summary(time.N)
  anova(time.N)
  time.N$coefficients
  check_model(time.N)
  
  ggplot(excr %>% filter(Species.code == 'YP'), 
         aes(x = Time.elapsed, y = log10(massnorm.P.excr))) +
    geom_point(size = 2) +
    geom_smooth(method = lm) +
    theme_classic()
  
  time.P <- lm(log10(massnorm.P.excr) ~ Time.elapsed, data = excr)
  summary(time.P)
  anova(time.P)
  time.P$coefficients
  check_model(time.P)
  
  ggplot(excr, aes(x = Time.elapsed, y = Temp)) +
    geom_point(size = 2) +
    #geom_smooth(method = lm) +
    theme_classic()
  
  # What about fish mass and incubation time or ambient DOC?
  
  ggplot(excr, aes(x = Time.elapsed, y = log10(Mass))) +
    geom_point(size = 2) +
    #geom_smooth(method = lm) +
    theme_classic()
  
  ggplot(excr, aes(x = AmDOC, y = log10(Mass))) +
    geom_point(size = 2) +
    geom_smooth() +
    theme_classic()
  
  ggplot(excr, aes(x = AmDOC, y = Temp)) +
    geom_point(size = 2) +
    geom_smooth() +
    theme_classic()
  
  
  


  
  
  
  