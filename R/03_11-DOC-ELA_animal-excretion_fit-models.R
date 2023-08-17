  #### fish nutrient and DOC excretion across a lake DOC gradient ####
  # This code was created by S. Klemet-N'Guessan in 2022 and 2023
  # R version 4.3.0
 
  ############################ MODELS ############################################# 
  # load libraries ----
  library(mgcv) # for GAM analysis (gam function)
  library(gratia) # to draw GAM
  library(vegan) # for NMDS
  library(car) # for Anova function
  library(rstatix) # for many pipe-friendly stat tools
  
  # PCA ----
  excr.pca <- excr.pca %>% select(-c(AmC1:AmC7))
  # change variable names for plotting purposes
  excr.pca.p <- excr.pca %>% rename_DOM()
  # do PCA
  pca <- princomp(excr.pca.p[, 3:13], cor = TRUE, scores = TRUE)
  biplot(pca)
  # check loadings
  pca$loadings
  # Extract PC1 scores
  pca.PC1 <- pca$scores[,1]
  # add it to current datasets
  excr.pca <- excr.pca %>% mutate(PC1 = pca.PC1)
  excr <- left_join(excr, excr.pca)
  
  # GAM ----
  # without trophic position and with population averages
  # create functions to calculate all HGAMs following a template
  # option 1
  hgam_s <- function(y, x, k1, df) {
    mod <- gam(y ~ s(x, k = k1, m = 2, bs = 'tp'),
               method = 'REML', data = df,
               family = tw())
    return(mod)
  }
  
  hgam <- function(y, x, k1, k2) {
    mod <- gam(y ~ s(x, by = Trophic.position,
                     k = k1, m = 2, bs = 'tp') +
                 s(Trophic.position, bs = 're', k = k2),
               method = 'REML', data = excr,
               family = tw())
    return(mod)
  }
  
  hgam_log <- function(y, x, k1, k2) {
    mod <- gam(y ~ s(x, by = Trophic.position,
                     k = k1, m = 2, bs = 'tp') +
                 s(Trophic.position, bs = 're', k = k2),
               method = 'REML', data = excr)
    return(mod)
  }
  
  hgam_null <- function(y, Trophic.position, k) {
    mod <- gam(y ~ s(Trophic.position, bs = 're', k = k), 
               method = 'REML', select = T, data = excr,
               family  = tw())
    return(mod)
  }
  
  hgam_null_log <- function(y, Trophic.position, k) {
    mod <- gam(y ~ s(Trophic.position, bs = 're', k = k), 
               method = 'REML', select = T, data = excr)
    return(mod)
  }
  
  # list to check all model details using lapply()
  gam.details <- list(summary = summary.gam, 
                      gam.check = gam.check, 
                      appraise = appraise, 
                      draw = draw)
  # ...N excretion ----
  # DOC
  gamNDOC_s <- hgam_s(excr$massnorm.N.excr, excr$AmDOC, 9, excr)
  gamNDOC <- hgam(excr$massnorm.N.excr, excr$AmDOC, 6, 2)
  lapply(gam.details, function(f) f(gamNDOC))
  AIC(gamNDOC_s, gamNDOC)
  # DOM
  gamNDOM <- hgam(excr$massnorm.N.excr, excr$PC1, 7, 2)
  lapply(gam.details, function(f) f(gamNDOM))
  
  # ...P excretion ----
  # DOC
  gamPDOC_s <- hgam_s(excr$massnorm.P.excr, excr$AmDOC, 7, excr)
  gamPDOC <- hgam(excr$massnorm.P.excr, excr$AmDOC, 5, 2)
  lapply(gam.details, function(f) f(gamPDOC))
  AIC(gamNDOC_s, gamNDOC)
  # DOM
  gamPDOM <- hgam(excr$massnorm.P.excr, excr$PC1, 6, 2)
  lapply(gam.details, function(f) f(gamPDOM))
  
  # ...N:P excretion ----
  # DOC
  gamNPDOC <- hgam_log(excr$massnorm.NP.excr, excr$AmDOC, 5, 2)
  lapply(gam.details, function(f) f(gamNPDOC))
  # DOM
  gamNPDOM <- hgam_log(excr$massnorm.NP.excr, excr$PC1, 6, 2)
  lapply(gam.details, function(f) f(gamNPDOM))
  
  # create list of gam functions to store
  # gam_DOC_list <- list(c(gamNDOC, gamPDOC))
  # gam_DOM_list <- list(c(gamNDOM, gamPDOM))
  # gam_DOC_log_list <- list(c(gamNPDOC))
  # gam_DOM_log_list <- list(c(gamNPDOM))
  
  # ANOVA ----
  aov_C <- function(x) {
    aov(log10(x) ~ DOC.level, 
        data = excr)
  }
  aovC <- aov_C(excr$massnorm.C.excr)
  summary(aovC)
  
  aovCN <- aov_C(excr$massnorm.CN.excr)
  summary(aovCN)
  posthoc <- tukey_hsd(aovCN, conf.level = .95)

  aovCP <- aov_C(excr$massnorm.CP.excr)
  summary(aovCP)
  
  # t-test ----
  # NEED TO LOOP IT SOMEHOW
  excr.DOM %>% 
    filter(type == 'C4') %>% 
    t_test(massnorm.excr ~ 1, mu = 0,  alternative = "greater")
  
  # NMDS ----
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
    as_tibble(rownames = 'sample') %>% 
    dplyr::mutate(
      ID = excr.nmds$ID,
      Site.name = excr.nmds$Site.name,
      Species.code = excr.nmds$Species.code
    )
  
