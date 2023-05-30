  #### fish excretion across a DOC gradient ####
  # This code was created by S. Klemet-N'Guessan in 2022
  # R version 4.1.3
  
  # load libraries and read datasets ----
  
  library(tidyverse)
  library(RColorBrewer)
  library(ggpubr) # for ggarrange function that allows to put all graphs on 1 page
  library(mgcv) # for GAM analysis (gam function)
  library(gratia) # to draw GAM
  
  # upload datasets ----
  er <- read.csv('data/2022-07-11_11-DOC-lakes_Mastersheet.csv',
                 stringsAsFactors = F, na.strings = c("", "NA", "."), 
                 strip.white = TRUE, sep = ",")
  
  str(er) 
  head(er)
  unique(er$P.excretion.rate..ug.h.ind.)
  
  # clean data, rename and add variables ----
  
  excr <- er %>%
    rename(Mass = Dry.mass..g.,
           AmHisDOC = Ambient.historical.DOC..mg.C.L.,
           AmDOC = DOC..mg.C.L.,
           P.excretion.rate = P.excretion.rate..ug.h.ind.,
           N.excretion.rate = N.excretion.rate..ug.h.ind.) %>% 
    # filter(Comments != 'dead during experiment') %>%
    mutate(Trophic.position = ifelse(Species.code == 'YP', 'C',
                                     ifelse(Species.code == 'PD',
                                            'I', 'O')),
           Log10.mass = log10(Mass),
           Log10.P.excretion.rate = round(log10(P.excretion.rate), 6),
           Log10.N.excretion.rate = round(log10(N.excretion.rate), 6)) %>%
    filter(!is.na(Log10.N.excretion.rate),
           !is.na(Log10.P.excretion.rate))
  
  # look at mass vs. N/P excretion ----
  # Fathead minnows only because coeffs are
  # too high for N excr (>1) + too low for P excr (<0.3)
  # N excretion
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
  ggplot(excr %>%  filter(Species.code == 'FM',
                          Log10.P.excretion.rate > 0.5),
         aes(x = Log10.mass, y = Log10.P.excretion.rate)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    theme(legend.position = 'none')
  
  # all species
  # N excretion
  ggplot(excr,
         aes(x = Log10.mass, y = Log10.N.excretion.rate)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    theme_classic() #+
    #stat_cor()
      
  modelN <- lm(Log10.N.excretion.rate ~ Log10.mass, 
               data = excr) 
  Ncoeff.allsp <- modelN$coefficients["Log10.mass"]
  
  ggplot(excr,
         aes(x = Log10.mass, y = Log10.P.excretion.rate)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    theme(legend.position = 'none') +
    theme_classic() #+
    #stat_cor()
  
  modelP <- lm(Log10.P.excretion.rate ~ Log10.mass, 
               data = excr) 
  Pcoeff.allsp <- modelP$coefficients["Log10.mass"]
  
  # ..Sort observations by species ----
  obs.spsummary <- excr %>% 
    group_by(Species.code) %>% 
    # Count number of observations by group, put the count in a new column, "n" (that's what tally does")
    tally() %>% 
    # Now arrange by smallest number of observations to largest
    arrange(n)
  
  head(obs.spsummary)
  
  # let's take only species with 10 or more observations and carry on.
  # Here I create a new data.frame that only has species with 10 or more observations. 
  newdf.sp <- obs.spsummary %>% 
    filter(n > 11) %>% 
    left_join(excr, by = "Species.code") %>% 
    select(-10)
  
  # Now get unique species in this new df
  
  # ....scaling exponent b for each species for N/P excretions ----
  # what are the unique species
  species <- unique(newdf.sp$Species.code)
  nb.species <- length(species)
  
  results.spdf <- data.frame() # an empty dataframe to put results into
  
  for (i in 1:nb.species) {
    subdf <- newdf.sp %>% 
      filter(Species.code == species[i]) # equivalent to 'Species X'
    subdf
    
    subdf %>% 
      select(Log10.N.excretion.rate, 
             Log10.P.excretion.rate, Log10.mass) 
    
    
    modelN <- lm(Log10.N.excretion.rate~Log10.mass, data = subdf)
    modelP <- lm(Log10.P.excretion.rate~Log10.mass, data = subdf)
    result <- data.frame(Species.code = species[i],
                         b.coeff.N.excr.sp = modelN$coefficients["Log10.mass"],
                         b.coeff.P.excr.sp = modelP$coefficients["Log10.mass"],
                         stringsAsFactors = FALSE)
    results.spdf <- bind_rows(results.spdf, result)
    cat(species[i], '\n') # to check where are in loop
  }
  results.spdf # Here are all of your results in one data.frame. 
  excr <- left_join(excr, results.spdf, by = "Species.code")
  
  # ..do mass-corrected excretion rates calculations ----
  excr <- excr %>% 
    mutate(massnorm.N.excr = N.excretion.rate/(Mass^(Ncoeff.allsp)),
           massnorm.P.excr = P.excretion.rate/(Mass^(Pcoeff.allsp)),
           # massnorm.N.excr = N.excretion.rate/(Mass^(b.coeff.N.excr.sp)),
           # massnorm.P.excr = P.excretion.rate/(Mass^(b.coeff.P.excr.sp)),
           massnorm.NP.excr = (massnorm.N.excr/massnorm.P.excr)/(14/31))
  
  # ..make excr dataset with one entry for each excretion average ----
  excr.sp <- excr %>% group_by(Site.name, Species.code, Trophic.position) %>% 
    summarise(massnorm.N.excr.sp = mean(massnorm.N.excr, na.rm = T),
              massnorm.P.excr.sp = mean(massnorm.P.excr, na.rm = T),
              massnorm.NP.excr.sp = mean(massnorm.NP.excr, na.rm = T),
              Log10.massnorm.N.excr.sp = log10(massnorm.N.excr.sp),
              Log10.massnorm.P.excr.sp = log10(massnorm.P.excr.sp),
              Log10.massnorm.NP.excr.sp = log10(massnorm.NP.excr.sp),
              AmDOC = mean(AmDOC))
  
  excr.smry <- excr %>% filter(Species.code != 'CTL1',
                               Species.code != 'CTL2') %>% 
    group_by(Site.name) %>% 
    summarise(massnorm.N.excr.sp = mean(massnorm.N.excr, na.rm = T),
              massnorm.P.excr.sp = mean(massnorm.P.excr, na.rm = T),
              massnorm.NP.excr.sp = mean(massnorm.NP.excr, na.rm = T),
              Log10.massnorm.N.excr.sp = log10(massnorm.N.excr.sp),
              Log10.massnorm.P.excr.sp = log10(massnorm.P.excr.sp),
              Log10.massnorm.NP.excr.sp = log10(massnorm.NP.excr.sp),
              AmHisDOC = mean(AmHisDOC))
  
  #### PRELIMINARY RESULTS ####
  # labels
  Trophic.labels <- c('Invert/piscivore', 'Invertivore', 'Omnivore')
  Species.labels <- c('fathead minnow', 'pearl dace', 
                    'white sucker', 'yellow perch')
  
  # ..Trophic position ----
  # N excretion vs. Hist DOC relative to trophic position
  Nexcrtp.p <- ggplot(excr.sp %>%  filter(Species.code != 'CTL1',
                              Species.code != 'CTL2'), 
         aes(x = AmHisDOC, y = massnorm.N.excr.sp)) +
    geom_point(aes(color = Trophic.position), 
               size = 4) +
    theme_classic(base_size = 20) +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-corrected',
                             paste(N~excretion~(μg~N/g/h))))) +
    scale_color_brewer(palette = "Dark2",
                       name = 'Trophic position',
                       labels = Trophic.labels)
  Nexcrtp.p
  
  # P excretion vs. Hist DOC relative to trophic position
  Pexcrtp.p <- ggplot(excr.sp %>%  filter(Species.code != 'CTL1',
                             Species.code != 'CTL2'), 
         aes(x = AmHisDOC, y = massnorm.P.excr.sp)) +
    geom_point(aes(color = Trophic.position), 
               size = 4) +
    theme_classic(base_size = 20) +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-corrected',
                             paste(P~excretion~(μg~P/g/h))))) +
    scale_color_brewer(palette = "Dark2",
                       name = 'Trophic position',
                       labels = Trophic.labels)
  Pexcrtp.p
  
  # put them together
  ggarrange(Nexcrtp.p, Pexcrtp.p,
            ncol = 2, align = "h", common.legend = T)
  
  ggsave('preliminary results/N_Pexcretion_tp.png', 
         width = 14, height = 7, 
         units = 'in', dpi = 600)
  
  # ..Species ----
  # N excretion vs. Hist DOC relative to species
  Nexcrsp.p <- ggplot(excr.sp %>%  filter(Species.code != 'CTL1',
                                          Species.code != 'CTL2'), 
                      aes(x = AmDOC, y = massnorm.N.excr.sp)) +
    geom_point(aes(color = Species.code), 
               size = 4) +
    theme_classic(base_size = 20) +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-corrected',
                             paste(N~excretion~(μg~N/g/h))))) +
    scale_color_brewer(palette = "Set1",
                       name = 'Species',
                       labels = Species.labels)
  Nexcrsp.p
  
  # P excretion vs. Hist DOC
  Pexcrsp.p <- ggplot(excr.sp %>%  filter(Species.code != 'CTL1',
                                          Species.code != 'CTL2'), 
                      aes(x = AmDOC, y = massnorm.P.excr.sp)) +
    geom_point(aes(color = Species.code), 
               size = 4) +
    theme_classic(base_size = 20) +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-corrected',
                             paste(P~excretion~(μg~P/g/h))))) +
    scale_color_brewer(palette = "Set1",
                       name = 'Species',
                       labels = Species.labels)
  Pexcrsp.p
  
  # put them together
  ggarrange(Nexcrsp.p, Pexcrsp.p,
            ncol = 2, align = "h", common.legend = T)
  
  ggsave('preliminary results/N_Pexcretion_sp.png', 
         width = 14, height = 7, 
         units = 'in', dpi = 600)
  
  # ..fish average ----
  # N excretion vs. Hist DOC
  Nexcrav.p <- ggplot(excr.smry, 
         aes(x = AmHisDOC, y = massnorm.N.excr.sp)) +
    geom_point(size = 4) +
    theme_classic(base_size = 20) +
  labs(x = 'DOC (mg C/L)',
       y = expression(atop('Mass-corrected',
                           paste(N~excretion~(μg~N/g/h)))))
  Nexcrav.p
  
  # P excretion vs. Hist DOC
  Pexcrav.p <- ggplot(excr, 
                       aes(x = AmDOC, y = massnorm.P.excr)) +
    geom_point(size = 4) +
    theme_classic(base_size = 20) +
    labs(x = 'DOC (mg C/L)',
         y = expression(atop('Mass-corrected',
                             paste(P~excretion~(μg~P/g/h))))) 
  Pexcrav.p
  
  # put them together
  ggarrange(Nexcrav.p, Pexcrav.p,
            ncol = 2, align = "h", common.legend = T)
  
  ggsave('preliminary results/N_Pexcretion_av.png', 
         width = 14, height = 7, 
         units = 'in', dpi = 600)
  