  #### fish nutrient and DOC excretion across a lake DOC gradient ####
  # This code was created by S. Klemet-N'Guessan in 2022 and 2023
  # R version 4.3.0
 
  ############################ MODELS ############################################# 
  # load libraries ----
  library(mgcv) # for GAM analysis (gam function)
  library(gratia) # to draw GAM
  library(vegan) # for NMDS
  library(cluster) # loaded for pairwiseAdonis
  library(pairwiseAdonis) # for NMDS posthoc
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
  hgam_s <- function(y, x, k1) {
    mod <- gam(y ~ s(x, k = k1, m = 2, bs = 'tp'),
               method = 'REML', data = excr,
               family = tw())
    return(mod)
  }
  
  hgam <- function(y, x, k1) {
    mod <- gam(y ~ s(x, k = k1, m = 2, bs = 'tp'),
               method = 'REML', data = excr,
               family = tw())
    return(mod)
  }
  
  hgam_log <- function(y, x, k1) {
    mod <- gam(log10(y) ~ s(x, k = k1, m = 2, bs = 'tp'),
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
  # gamNDOC_s <- hgam_s(excr$massnorm.N.excr, excr$AmDOC, 9, excr)
  gamNDOC <- hgam(excr$massnorm.N.excr, excr$AmDOC, 9)
  lapply(gam.details, function(f) f(gamNDOC))
  #AIC(gamNDOC_s, gamNDOC)
  # DOM
  gamNDOM <- hgam(excr$massnorm.N.excr, excr$PC1, 9)
  lapply(gam.details, function(f) f(gamNDOM))
  
  # ...P excretion ----
  # DOC
  # gamPDOC_s <- hgam_s(excr$massnorm.P.excr, excr$AmDOC, 7, excr)
  gamPDOC <- hgam(excr$massnorm.P.excr, excr$AmDOC, 5)
  lapply(gam.details, function(f) f(gamPDOC))
  #AIC(gamNDOC_s, gamNDOC)
  # DOM
  gamPDOM <- hgam(excr$massnorm.P.excr, excr$PC1, 6)
  lapply(gam.details, function(f) f(gamPDOM))
  
  # ...N:P excretion ----
  # DOC
  gamNPDOC <- hgam_log(excr$massnorm.NP.excr, excr$AmDOC, 6)
  lapply(gam.details, function(f) f(gamNPDOC))
  # DOM
  gamNPDOM <- hgam_log(excr$massnorm.NP.excr, excr$PC1, 7)
  lapply(gam.details, function(f) f(gamNPDOM))
  
  # create list of gam functions to store
  # gam_DOC_list <- list(c(gamNDOC, gamPDOC))
  # gam_DOM_list <- list(c(gamNDOM, gamPDOM))
  # gam_DOC_log_list <- list(c(gamNPDOC))
  # gam_DOM_log_list <- list(c(gamNPDOM))
  
  # ANOVA ----
  excr.aov <- excr %>% filter(!is.na(massnorm.C.excr),
                              massnorm.C.excr < 32.56)
  aov_DOC <- function(x) {
    m <- aov(log10(x) ~ DOC.level, 
        data = excr.aov)
    print(summary(m))
    return(m)
  }
  
  aovC <- aov_DOC(excr.aov$massnorm.C.excr)
  aovC.posthoc <- tukey_hsd(aovC) 
  aovC.posthoc
  
  aovCN <- aov_DOC(excr.aov$massnorm.CN.excr)
  aovCN.posthoc <- tukey_hsd(aovCN) 
  aovCN.posthoc
  
  aovCP <- aov_DOC(excr.aov$massnorm.CP.excr)
  
  # t-test ----
  # NEED TO LOOP IT SOMEHOW
  excr.DOM %>% 
    filter(type == 'C4') %>% 
    t_test(massnorm.excr ~ 1, mu = 0,  alternative = "greater")
  
  # NMDS ----
  # separate dataset between each site (low, medium, high)
  excr.nmds.l <- excr.nmds %>% 
    dplyr::filter(Site.name == 'L224') 
  excr.nmds.m <- excr.nmds %>% 
    dplyr::filter(Site.name == 'L239') 
  excr.nmds.h <- excr.nmds %>% 
    dplyr::filter(Site.name == 'L222')
  
  # transform NMDS dataset to a matrix
  mtx <- function(df) {
    df %>% select(-c(ID, Source, Site.name, Trophic.position)) %>% 
      as.matrix()
  }
  excr.nmds.lm <- mtx(excr.nmds.l)
  excr.nmds.mm <- mtx(excr.nmds.m)
  excr.nmds.hm <- mtx(excr.nmds.h)
  excr.nmds.allm <- mtx(excr.nmds)
  
  # do NMDS
  nmds <- function(df) {
    set.seed(1)
    m <- metaMDS(df, autotransform = FALSE,
                 distance = "bray",
                 engine = "monoMDS",
                 k = 2,
                 weakties = TRUE,
                 model = "global",
                 maxit = 300,
                 try = 40,
                 trymax = 100)
    stressplot(m)
    cat('stress:', m$stress)
    return(m)
  }
  
  nmds.l <- nmds(excr.nmds.lm)
  nmds.m <- nmds(excr.nmds.mm)
  nmds.h <- nmds(excr.nmds.hm)
  nmds.all <- nmds(excr.nmds.allm)
  
  # First create a data frame of the scores from the individual sites.
  # This data frame will contain x and y values for where sites are located.
  nmds_score <- function(m, df) {
    scores(m)$sites %>% 
      as_tibble(rownames = 'sample') %>% 
      dplyr::mutate(
        ID = df$ID,
        Site.name = df$Site.name,
        Source = df$Source,
        Trophic.position = df$Trophic.position
      )
  }
  
  nmds.l.scores <- nmds_score(nmds.l, excr.nmds.l)
  nmds.m.scores <- nmds_score(nmds.m, excr.nmds.m)
  nmds.h.scores <- nmds_score(nmds.h, excr.nmds.h)
  nmds.all.scores <- nmds_score(nmds.all, excr.nmds)
  
  # # get hull data to draw polygons on plot
  # FM <-
  #   nmds.l.scores[nmds.l.scores$Source == "FM",][chull(nmds.l.scores[nmds.l.scores$Source ==
  #                                                                      "FM", c("NMDS1", "NMDS2")]),]
  # PD <-
  #   nmds.l.scores[nmds.l.scores$Source == "PD",][chull(nmds.l.scores[nmds.l.scores$Source ==
  #                                                                      "PD", c("NMDS1", "NMDS2")]),]
  # WS <-
  #   nmds.l.scores[nmds.l.scores$Source == "WS",][chull(nmds.l.scores[nmds.l.scores$Source ==
  #                                                                      "WS", c("NMDS1", "NMDS2")]),]
  # L224 <-
  #   nmds.l.scores[nmds.l.scores$Source == "AmL224",][chull(nmds.l.scores[nmds.l.scores$Source ==
  # "AmL224", c("NMDS1", "NMDS2")]),]
  
  O <-
    nmds.l.scores[nmds.l.scores$Trophic.position == "O",][chull(nmds.l.scores[nmds.l.scores$Trophic.position ==
                                                                       "O", c("NMDS1", "NMDS2")]),]
  I <-
    nmds.l.scores[nmds.l.scores$Trophic.position == "I",][chull(nmds.l.scores[nmds.l.scores$Trophic.position ==
                                                                       "I", c("NMDS1", "NMDS2")]),]
  L224 <-
    nmds.l.scores[nmds.l.scores$Trophic.position == "AmL224",][chull(nmds.l.scores[nmds.l.scores$Trophic.position ==
                                                                           "AmL224", c("NMDS1", "NMDS2")]),]
  # hull.l <- rbind(FM, PD, WS, L224)
  hull.l <- rbind(I, O, L224)
  
  # do PERMANOVA to test for differences between species/ambient water
  perma <- function(matrix, df) {
    set.seed(1)
    permanova <- adonis2(matrix ~ Trophic.position, df)
    print(permanova)
  }
  perma.l <- perma(excr.nmds.lm, excr.nmds.l)
  posthoc.l <- pairwise.adonis2(excr.nmds.lm ~ Trophic.position, excr.nmds.l)
  posthoc.l
  perma.m <- perma(excr.nmds.mm, excr.nmds.m)
  perma.h <- perma(excr.nmds.hm, excr.nmds.h)
  perma.all <- adonis2(excr.nmds.allm ~ Trophic.position*Site.name, excr.nmds, 
                       methods = 'bray')
  perma.all
  posthoc.all <- pairwise.adonis2(excr.nmds.allm ~ Trophic.position, excr.nmds)
  posthoc.all
  
  # do posthoc analysis
  
  # calculate homogeneity of group dispersion (variance)
  disper <- function(matrix, var){
    dis <- vegdist(matrix)
    result <- betadisper(dis, var)
    print(anova(result))
  }
  # there is heterogeneity of variances for all tests
  disper.l <- disper(excr.nmds.lm, excr.nmds.l$Source)
  disper.m <- disper(excr.nmds.mm, excr.nmds.m$Source)
  disper.h <- disper(excr.nmds.hm, excr.nmds.h$Source)
  disper.all <- disper(excr.nmds.allm, excr.nmds$Source)
  