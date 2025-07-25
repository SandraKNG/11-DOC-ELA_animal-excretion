  #### Fish supply distinct nutrients and dissolved organic matter composition ####
  ### relative to ambient concentrations in northern lakes ####
  
  # This code was created by S. Klemet-N'Guessan in 2022-2024
  # R version 4.3.0
 
  ############################ MODELS ############################################# 
  # load libraries ----
  library(PCAtest)
  library(mgcv) # for GAM analysis (gam function)
  library(gratia) # to draw GAM
  library(vegan) # for NMDS
  library(cluster) # loaded for pairwiseAdonis
  library(pairwiseAdonis) # for NMDS posthoc
  library(car) # for Anova function
  library(rstatix) # for many pipe-friendly stat tools
  library(writexl) # to export excel file
  library(ggcorrplot)
  library(PCAtest)
  
  # PCA ----
  # PCA to map all lake parameters except for DOM composition
  excr.pca.all <- excr.pca %>% 
    select(c(Site.name:AmTDN, Area:Part.P), 
           -c(pH, Zmax, Thermo.depth, Conductivity)) %>% 
    rename(DOC = AmDOC,
           TDP = AmTDP,
           TDN = AmTDN) 
  # do PCA
  pca.all <- princomp(excr.pca.all[, 2:11], cor = TRUE, scores = TRUE)
  biplot(pca.all)
  # check loadings
  pca.all$loadings
  pca.result <- PCAtest(excr.pca.all[, 2:11], 100, 100, 0.05, varcorr = FALSE, 
                    counter = FALSE, plot = TRUE)
  # Extract PC1 scores
  pca.PC1 <- pca.all$scores[,1]
  # add it to current datasets
  excr.pca <- excr.pca %>% mutate(PC1_all = pca.PC1)
  excr <- left_join(excr, excr.pca)
  
  # PCA to extract PC1 for DOM composition
  # change variable names for plotting purposes
  excr.pca.DOM <- excr.pca %>% 
    select(-c(AmC1:Part.P, C_humicper:PC1_all)) %>% 
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
  cor.test(excr.pca$AmDOC, excr.pca$AmTDP)
  cor.test(excr.pca$AmDOC, excr.pca$AmTDN)
  cor.test(excr.pca$AmDOC, excr.pca$Zmean)
  
  # correlation matrix
  excr.cor <- excr.pca %>% select(-c(AmC1:AmC7))
  cor.mx <- round(cor(excr.cor[2:26], use = "pairwise.complete.obs"), 1)
  # Computing correlation matrix with p-values 
  corrp.mat <- cor_pmat(excr.cor[2:26], use = "pairwise.complete.obs") 
  
  
  # GAM ----
  # without trophic position and with population averages
  # create functions to calculate all HGAMs following a template
  
  hgam <- function(y, k) {
    mod <- gam(y ~ s(AmDOC, by = Trophic.position2, k = k, m = 2, bs = 'tp') +
                 s(Trophic.position2, bs = 're'),
               method = 'REML', data = excr,
               family = tw())
    return(mod)
  }
#
  hgam_log <- function(y, k) {
    mod <- gam(log10(y) ~ s(AmDOC, by = Trophic.position2, k = k, m = 2, bs = 'tp') +
                 s(Trophic.position2, bs = 're'), 
               method = 'REML', data = excr)
    return(mod)
  }
  
  hgam_null <- function(y) {
    mod <- gam(y ~ s(Trophic.position2, bs = 're'),
               method = 'REML',  data = excr,
               family  = tw())
    return(mod)
  }
  
  hgam_null_log <- function(y) {
    mod <- gam(log10(y) ~ s(Trophic.position2, bs = 're'),
               method = 'REML', data = excr)
    return(mod)
  }

  hgam_lake <- function(y) {
    mod <- gam(y ~ Site.name +
                 s(Trophic.position2, bs = 're'),
               method = 'REML', data = excr,
               family  = tw())
    return(mod)
  }
  
  hgam_lake_log <- function(y) {
    mod <- gam(log10(y) ~ Site.name +
                 s(Trophic.position2, bs = 're'),
               method = 'REML', data = excr)
    return(mod)
  }
  
  # list to check all model details using lapply()
  gam.details <- list(summary = summary.gam, 
                      gam.check = gam.check, 
                      appraise = appraise, 
                      draw = draw)
  
  # ...N excretion ----
  # DOC
  gamNDOC <- hgam(excr$massnorm.N.excr, 6)
  lapply(gam.details, function(f) f(gamNDOC))
  gamN.null <- hgam_null(excr$massnorm.N.excr)
  gamN.lake <- hgam_lake(excr$massnorm.N.excr)
  summary(gamN.lake)
  AIC(gamNDOC)
  AIC(gamN.null)
   
  # ...P excretion ----
  # DOC
  gamPDOC <- hgam(excr$massnorm.P.excr, 5)
  lapply(gam.details, function(f) f(gamPDOC))
  gamP.null <- hgam_null(excr$massnorm.P.excr)
  gamP.lake <- hgam_lake(excr$massnorm.P.excr)
  AIC(gamPDOC)
  
  # ...N:P excretion ----
  # DOC
  gamNPDOC <- hgam_log(excr$massnorm.NP.excr, 5)
  lapply(gam.details, function(f) f(gamNPDOC))
  gamNP.null <- hgam_null_log(excr$massnorm.NP.excr)
  gamNP.lake <- hgam_lake_log(excr$massnorm.NP.excr)
  AIC(gamNPDOC)
  
  # make AIC table ----
  AIC.tbl <- AIC(gamNDOC, gamN.lake, gamN.null,  
                 gamPDOC, gamP.lake, gamP.null,
                 gamNPDOC, gamNP.lake, gamNP.null)
  
  AIC.table <- AIC.tbl %>%  rownames_to_column(var = "Model") %>%
    mutate(df = round(df, digits = 0),
           AIC = round(AIC, digits = 0))
  write_xlsx(AIC.table, 'output/AIC_table.xlsx')
  
  # ANOVA ----
  excr.aov <- excr %>%
    filter(!is.na(massnorm.C.excr),
                              massnorm.C.excr < 32.56) %>% 
    mutate(across(where(is.numeric),
           ~ if_else(. < 0, NA, .)))
  
  aov_DOCM <- function(x) {
    m <- lm(log10(x) ~ DOC.level, data = excr.aov)
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
  anova.results <- bind_rows(anova_results, .id = "response")
  tukey.results <- bind_rows(tukey_results, .id = "response")
  tukey.results <- tukey.results %>% 
    arrange(p.adj, group1) 
  
    
  print(tukey.results, n = 33)
  
  # t-test ----
  excr.ttest <- excr.var 
  t_test_results <- list()
  
  # Loop through selected columns
  for (col in selected_cols) {
    result <- excr.ttest %>%
      t_test(as.formula(paste(col, "~ 1")), mu = 0, alternative = "greater")
    
    # Store the result in the list
    t_test_results[[col]] <- result
    
    # Print progress
    cat("Processed column:", col, "\n")
  }
  
  # Convert the list of results to a tibble
  t.test.results <- bind_rows(t_test_results) 
  t.test.results <- t.test.results %>% 
    arrange(p) %>% 
    rename(response = .y.)
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
    df %>% select(-c(ID, Source, Site.name, Trophic.position2)) %>% 
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
        Trophic.position2 = df$Trophic.position2
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
    permanova <- adonis2(matrix ~ Trophic.position2, df)
    print(permanova)
  }
  
  perma.l <- perma(excr.nmds.lm, excr.nmds.l)
  posthoc.l <- pairwise.adonis2(excr.nmds.lm ~ Trophic.position2, excr.nmds.l)
  posthoc.l
  perma.m <- perma(excr.nmds.mm, excr.nmds.m)
  perma.h <- perma(excr.nmds.hm, excr.nmds.h)
  perma.all <- adonis2(excr.nmds.allm ~ Source, excr.nmds, 
                       methods = 'bray')
  perma.all
  posthoc.all <- pairwise.adonis2(excr.nmds.allm ~ Source, excr.nmds)
  posthoc.all
  
  # calculate homogeneity of group dispersion (variance)
  disper <- function(matrix, var){
    dis <- vegdist(matrix)
    result <- betadisper(dis, var)
    print(anova(result))
  }
  
  # there is heterogeneity of variances for all tests
  disper.l <- disper(excr.nmds.lm, excr.nmds.l$Trophic.position2)
  disper.m <- disper(excr.nmds.mm, excr.nmds.m$Trophic.position2)
  disper.h <- disper(excr.nmds.hm, excr.nmds.h$Trophic.position2)
  disper.all <- disper(excr.nmds.allm, excr.nmds$Source)

  # export test result tables ----
  write_csv(pca.result, "output/PCAtest_results.csv")
  write_csv(tukey.results, "output/DOMexcr_tukey_results.csv")
  write_csv(t.test.results, "output/DOMexcr_ttest_results.csv")
  