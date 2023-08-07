  #### *fish nutrient and DOC excretion across a lake DOC gradient ####
  # This code was created by S. Klemet-N'Guessan in 2022 and 2023
  # R version 4.3.0
  
  # load libraries and read datasets ----
  
  library(tidyverse)
  library(RColorBrewer)
  library(MetBrewer) # Met brewer color palette
  library(datawizard) # to do summary statistics
  
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
           AmDOC = Ambient.DOC..mg.C.L.,
           P.excretion.rate = P.excretion.rate..ug.h.ind.,
           N.excretion.rate = N.excretion.rate..ug.h.ind.,
           C.excretion.rate = C.excretion..mg.h.ind.,
           AmSUVA = Ambient.SUVA,
           AmBA = Ambient.BA,
           AmSR = Ambient.SR,
           AmFI = Ambient.FI,
           AmHIX = Ambient.HIX.ohno,
           HIX.excretion = HIX.ohno.excretion) %>% 
    # filter(Comments != 'dead during experiment') %>%
    mutate(Trophic.position = factor(ifelse(Species.code == 'YP', 'C',
                                     ifelse(Species.code == 'PD',
                                            'I', 'O'))),
           DOC.level = factor(ifelse(AmDOC == 3.462, 1, 
                                     ifelse(AmDOC == 7.156, 2,
                                            ifelse(AmDOC == 10.04, 3, NA)))),
           Log10.mass = log10(Mass),
           Log10.P.excretion.rate = round(log10(P.excretion.rate), 6),
           Log10.N.excretion.rate = round(log10(N.excretion.rate), 6),
           Log10.C.excretion.rate = round(log10(C.excretion.rate), 6)) %>%
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
  
  # P excretion
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
  
  # C excretion
  ggplot(excr,
         aes(x = Log10.mass, y = Log10.C.excretion.rate)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    theme(legend.position = 'none') +
    theme_classic() #+
  #stat_cor()
  
  modelC <- lm(Log10.C.excretion.rate ~ Log10.mass, 
               data = excr) 
  Ccoeff.allsp <- modelC$coefficients["Log10.mass"]
  
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
           massnorm.C.excr = C.excretion.rate/(Mass^(Ccoeff.allsp)),
           # massnorm.N.excr = N.excretion.rate/(Mass^(b.coeff.N.excr.sp)),
           # massnorm.P.excr = P.excretion.rate/(Mass^(b.coeff.P.excr.sp)),
           massnorm.NP.excr = (massnorm.N.excr/massnorm.P.excr)/(14/31),
           massnorm.CN.excr = (massnorm.C.excr/massnorm.N.excr)/(12/14),
           massnorm.CP.excr = (massnorm.C.excr/massnorm.P.excr)/(12/31))
  
  # ..make excr dataset with one entry for each excretion average ----
  # N + P excretion species average
  excr.sp <- excr %>% group_by(Site.name, Species.code, Trophic.position) %>% 
    reframe(massnorm.N.excr.sp = mean(massnorm.N.excr, na.rm = T),
            massnorm.P.excr.sp = mean(massnorm.P.excr, na.rm = T),
            massnorm.NP.excr.sp = mean(massnorm.NP.excr, na.rm = T),
            Log10.massnorm.N.excr.sp = log10(massnorm.N.excr.sp),
            Log10.massnorm.P.excr.sp = log10(massnorm.P.excr.sp),
            Log10.massnorm.NP.excr.sp = log10(massnorm.NP.excr.sp),
            AmDOC = mean(AmDOC),
            AmBA = mean(AmBA))
  
  # C excretion species average
  excr.C.sp <- excr %>% filter(C.excretion.rate > 0) %>% 
    group_by(Site.name, Species.code, Trophic.position ) %>% 
    reframe(massnorm.C.excr.sp = mean(massnorm.C.excr, na.rm = T),
      massnorm.NP.excr.sp = mean(massnorm.NP.excr, na.rm = T),
      massnorm.CN.excr.sp = mean(massnorm.CN.excr, na.rm = T),
      massnorm.CP.excr.sp = mean(massnorm.CP.excr, na.rm = T),
      AmDOC = mean(AmDOC)) 
  excr.C.sp <- excr.C.sp %>% 
    mutate(DOC.level = factor(ifelse(AmDOC == 3.462, 'low', 
                              ifelse(AmDOC == 7.156, 'med',
                                     ifelse(AmDOC == 10.04, 'high', NA)))))
  
  excr.smry <- excr %>% filter(Species.code != 'CTL1',
                               Species.code != 'CTL2',
                               Species.code != 'CTL3',
                               Species.code != 'CTL4') %>% 
    group_by(Site.name) %>% 
    reframe(massnorm.N.excr.sp = mean(massnorm.N.excr, na.rm = T),
              massnorm.P.excr.sp = mean(massnorm.P.excr, na.rm = T),
              massnorm.C.excr.sp = mean(massnorm.C.excr, na.rm = T),
              massnorm.NP.excr.sp = mean(massnorm.NP.excr, na.rm = T),
              massnorm.CN.excr.sp = mean(massnorm.CN.excr, na.rm = T),
              massnorm.CP.excr.sp = mean(massnorm.CP.excr, na.rm = T),
              Log10.massnorm.N.excr.sp = log10(massnorm.N.excr.sp),
              Log10.massnorm.P.excr.sp = log10(massnorm.P.excr.sp),
              Log10.massnorm.NP.excr.sp = log10(massnorm.NP.excr.sp),
              AmDOC = mean(AmDOC))

  # ..summary statistics ----
  excr.ss <- excr %>% 
    select(c('massnorm.N.excr', 'massnorm.P.excr', 'massnorm.NP.excr',
             'Mass')) %>% 
    describe_distribution()
  
  excrtp.ss <- excr %>% group_by(Trophic.position) %>%
    select(c('massnorm.N.excr', 'massnorm.P.excr', 'massnorm.NP.excr',
             'Mass')) %>% 
    describe_distribution()