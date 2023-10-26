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
  # PCA to map all lake parameters except for DOM composition
  excr.pca.all <- excr.pca %>% 
    select(c(Site.name:AmTDP, Area:Part.P), 
           -c(pH, Zmean, Zmax)) %>% 
    rename(DOC = AmDOC,
           TDP = AmTDP) %>% 
    mutate(Thermo.depth = replace_na(Thermo.depth, 1.70))
  # do PCA
  pca.all <- princomp(excr.pca.all[, 2:11],cor = TRUE, scores = TRUE)
  biplot(pca.all)
  # check loadings
  pca.all$loadings
  # Extract PC1 scores
  pca.PC1 <- pca.all$scores[,1]
  # add it to current datasets
  excr.pca <- excr.pca %>% mutate(PC1_all = pca.PC1)
  excr <- left_join(excr, excr.pca)
  
  # PCA to extract PC1 for DOM composition
  # change variable names for plotting purposes
  excr.pca.DOM <- excr.pca %>% 
    select(-c(AmC1:Part.P)) %>% 
    rename_DOM() %>% 
    rename_with(~ gsub("per", "", .x, fixed = TRUE)) 
  # do PCA
  pca.DOM <- princomp(excr.pca.DOM[, 5:15], cor = TRUE, scores = TRUE)
  biplot(pca.DOM)
  # check loadings
  pca.DOM$loadings
  # Extract PC1 scores
  pca.DOM.PC1 <- pca.DOM$scores[,1]
  # add it to current datasets
  excr.pca <- excr.pca %>% mutate(PC1 = pca.DOM.PC1)
  excr <- left_join(excr, excr.pca)
  
  # correlation tests ----
  cor.test(excr.pca$AmDOC, excr.pca$PC1)
  ggplot(excr.pca, aes(AmDOC, PC1)) +
    geom_point(aes(colour = Site.name), size = 4) +
    geom_smooth(method = lm, colour = 'black') +
    xlab('DOC (mg C/L)') +
    theme_classic(base_size = 20) +
    annotate('text', x = 4.5, y = 7, 
             size = 5, label = 'cor = 0.9, p < 0.001')
  cor.test(excr.pca$AmDOC, excr.pca$AmTDP)
  
  # correlation matrix
  excr.cor <- excr.pca %>% select(-c(AmC1:AmC7))
  cor.mx <- cor(excr.cor[2:26])
  corrplot(cor.mx)
  
  # GAM ----
  # without trophic position and with population averages
  # create functions to calculate all HGAMs following a template
  
  hgam <- function(y, k1) {
    mod <- gam(y ~ s(AmDOC, by = Trophic.position2, k = k1, m = 2, bs = 'tp') +
                 s(Trophic.position2, bs = 're', k = 3),
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
  gamNDOC_s <- hgam(excr$massnorm.N.excr, 9)
  lapply(gam.details, function(f) f(gamNDOC_s))
  AIC(gamNDOC_s, gamNDOC)
  # # DOM
  # gamNDOM <- hgam(excr$massnorm.N.excr, excr$PC1, 9)
  # lapply(gam.details, function(f) f(gamNDOM))
  
  # ...P excretion ----
  # DOC
  # gamPDOC_s <- hgam_s(excr$massnorm.P.excr, excr$AmDOC, 7, excr)
  gamPDOC <- hgam(excr$massnorm.P.excr, 5)
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
  
  # ANOVA ----
  excr.aov <- excr %>% filter(!is.na(massnorm.C.excr),
                              massnorm.C.excr < 32.56)
  aov_DOCM <- function(x) {
    m <- aov(log10(x)~ DOC.level, data = excr.aov)
    print(summary(m))
    return(m)
  }
  
  # DOC tests
  aovC <- aov_DOCM(excr.aov$massnorm.C.excr)
  aovC.posthoc <- tukey_hsd(aovC) 
  aovC.posthoc
  
  aovCN <- aov_DOCM(excr.aov$massnorm.CN.excr)
  aovCN.posthoc <- tukey_hsd(aovCN) 
  aovCN.posthoc
  
  aovCP <- aov_DOCM(excr.aov$massnorm.CP.excr)
  
  # DOM tests
  # Get column indices between "massnorm.SUVA.excr" and "massnorm.C7.excr"
  start_col <- which(names(excr.var) == "massnorm.SUVA.excr")
  end_col <- which(names(excr.var) == "massnorm.C7.excr")
  selected_cols <- names(excr.var)[(start_col):(end_col)]
  
  # Create an empty list to store aov results
  anova_results <- list()
  tukey_results <- list()
  
  for (col in selected_cols) {
    aov_result <- Anova(aov_DOCM(excr.aov[[col]]))
    tukey_result <- tukey_hsd(excr.aov, log10(excr.aov[[col]]) ~ DOC.level)
  
  # Store the result in the list
  anova_results[[col]] <- aov_result
  tukey_results[[col]] <- tukey_result
  
  # Print progress
  cat("Processed column:", col, "\n")
  }
  
  # Convert the list of results to tibbles
  anova.results <- bind_rows(anova_results, .id = "column_name")
  tukey.results <- bind_rows(tukey_results, .id = "column_name")
  tukey.results <- tukey.results %>% arrange(p.adj.signif, column_name)
  print(tukey.results, n = 33)
  
  # t-test ----
  # Create an empty list to store t-test results
  t_test_results <- list()
  
  # Loop through selected columns
  for (col in selected_cols) {
    result <- excr.var %>%
      t_test(as.formula(paste(col, "~ 1")), mu = 0, alternative = "greater")
    
    # Store the result in the list
    t_test_results[[col]] <- result
    
    # Print progress
    cat("Processed column:", col, "\n")
  }
  
  # Convert the list of results to a tibble
  t.test.results <- bind_rows(t_test_results) 
  t.test.results <- t.test.results %>% arrange(p)
  t.test.results
  
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
  
  # do PERMANOVA to test for differences between species/ambient water
  # AND do posthoc analysis
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
  
  # calculate homogeneity of group dispersion (variance)
  disper <- function(matrix, var){
    dis <- vegdist(matrix)
    result <- betadisper(dis, var)
    print(anova(result))
  }
  
  # there is heterogeneity of variances for all tests
  disper.l <- disper(excr.nmds.lm, excr.nmds.l$Trophic.position)
  disper.m <- disper(excr.nmds.mm, excr.nmds.m$Trophic.position)
  disper.h <- disper(excr.nmds.hm, excr.nmds.h$Trophic.position)
  disper.all <- disper(excr.nmds.allm, excr.nmds$Trophic.position)
  