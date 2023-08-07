  #### fish nutrient and DOC excretion across a lake DOC gradient ####
  # This code was created by S. Klemet-N'Guessan in 2022 and 2023
  # R version 4.3.0
 
  ############################ MODELS ############################################# 
  # load libraries ----
  library(mgcv) # for GAM analysis (gam function)
  library(gratia) # to draw GAM
  library(vegan) # for NMDS
  library(car) # for Anova function
  
  # do GAM ----
  # without trophic position and with population averages
  gamNDOC <- gam(massnorm.N.excr.sp ~ s(AmDOC,
                                        k = 7, m = 2, bs = 'tp'),
                    method = 'REML', data = excr.sp)
  summary(gamNDOC)
  gam.check(gamNDOC)
  appraise(gamNDOC, method = 'simulate')
  draw(gamNDOC)
  
  # P excretion
  gamPDOC <- gam(massnorm.P.excr.sp ~ s(AmDOC, 
                                        k = 5, m = 2, bs = 'tp'),
                    method = 'REML', data = excr.sp)
  summary(gamPDOC)
  gam.check(gamPDOC)
  appraise(gamPDOC, method = 'simulate')
  draw(gamPDOC)
  
  # N:P excretion
  gamNPDOC <- gam(log10(massnorm.NP.excr.sp) ~ s(AmDOC, 
                                        k = 5, m = 2, bs = 'tp'),
                 method = 'REML', data = excr.sp)
  summary(gamNPDOC)
  gam.check(gamNPDOC)
  appraise(gamNPDOC, method = 'simulate')
  draw(gamNPDOC)
  
  # with trophic position
  gamNDOC.tp <- gam(massnorm.N.excr ~ s(AmDOC, by = Trophic.position,
                                     k = 7, m = 2, bs = 'tp') +
                   s(Trophic.position, bs = 're', k = 2),
                 method = 'REML', data = excr,
                 family = tw())
  summary(gamNDOC.tp)
  gam.check(gamNDOC.tp)
  appraise(gamNDOC.tp, method = 'simulate')
  draw(gamNDOC.tp)
  
  gamPDOC.tp <- gam(massnorm.P.excr ~ s(AmDOC, by = Trophic.position,
                                     k = 7, m = 2, bs = 'tp') +
                   s(Trophic.position, bs = 're', k = 2),
                 method = 'REML', data = excr,
                 family = tw())
  summary(gamPDOC.tp)
  gam.check(gamPDOC.tp)
  appraise(gamPDOC.tp, method = 'simulate')
  draw(gamPDOC.tp)
  
  # do LM ----
  lmN <- lm(massnorm.N.excr.sp ~ AmDOC, data = excr.sp)
  summary(lmN)
  Anova(lmN)
  
  lmP <- lm(massnorm.P.excr.sp ~ AmDOC, data = excr.sp)
  summary(lmP)
  Anova(lmP)
  
  # do ANOVA
  aov.C <- aov(log10(massnorm.C.excr) ~ DOC.level, 
               data = excr %>% filter(!is.na(DOC.level),
                                      C.excretion.rate > 0))
  Anova(aov.C)
  
  aov.CN <- aov(log10(massnorm.CN.excr) ~ DOC.level, 
                data = excr %>% filter(!is.na(DOC.level),
                                       C.excretion.rate > 0))
  Anova(aov.CN)
  
  aov.CP <- aov(log10(massnorm.CP.excr) ~ DOC.level, 
                data = excr %>% filter(!is.na(DOC.level),
                                       C.excretion.rate > 0))
  Anova(aov.CP)
