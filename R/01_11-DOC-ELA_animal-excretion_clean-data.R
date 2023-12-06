  #### ***fish nutrient and DOC excretion across a lake DOC gradient ####
  # This code was created by S. Klemet-N'Guessan in 2022 and 2023
  # R version 4.3.0
  
  # Load libraries and read datasets ----
  
  library(tidyverse)
  library(datawizard) # to do summary statistics
  
  # Upload datasets ----
  er <- read_csv('data/2022-07-11_11-DOC-lakes_Mastersheet.csv', 
                 na = c("NA", ""))
  
  er
  str(er)
  head(er)
  
  # Clean data, rename and add variables ----
  unique(er$Comments)
  excr <- er %>%
    rename(
      Mass = Dry.mass,
      AmHisDOC = Amhistorical.DOC,
      Ambient.HIX = AmHIX,
      AmSUVA = AmSUVA254,
      AmHIX = AmHIX.ohno
    ) %>% 
    mutate(
      Species.code = factor(Species.code),
      Trophic.position = factor(
        ifelse(Species.code == 'YP', 'C', ifelse(Species.code == 'PD',
                                                 'I', 'O')
      )),
      Trophic.position2 = factor(
        ifelse(Species.code == 'YP' | Species.code == 'PD', 'C', 'O')
      )
      ) %>% 
    # filter out control samples + fish that could have had a contaminated water
    filter(
      !Species.code %in% c('CTL1', 'CTL2', 'CTL3', 'CTL4'),
      !Comments %in% c('dead during experiment', 
                      'water reddish when at filtering stage',
                      "put net from anesthetic bin in bag")) %>% 
    mutate(
      DOC.level = factor(ifelse(Site.name == 'L224', 'low', 
                                ifelse(Site.name == 'L239', 'med',
                                       ifelse(Site.name == 'L222', 'high', NA)))),
      DOC.level = fct_relevel(DOC.level, c('low', 'med', 'high'))
    )

  # Look at N/P excretion vs. mass ----
  # ...Fathead minnows only  ----
  # because coeffs are too high for N excr (>1) + too low for P excr (<0.3)
  # N excretion
  # ggplot(excr %>%  filter(Species.code == 'FM',
  #                         Log10.N.excretion.rate > 0.5),
  #        aes(x = Log10.mass, y = Log10.N.excretion.rate)) +
  #   geom_point() +
  #   geom_smooth(method = 'lm') +
  #   theme(legend.position = 'none')
  # modelN.FM <- lm(Log10.N.excretion.rate ~ Log10.mass, 
  #              data = excr %>%  
  #                filter(Species.code == 'FM',
  #                       Log10.N.excretion.rate > 0.5))
  # modelN.FM$coefficients["Log10.mass"]
  # 
  # # P excretion
  # ggplot(excr %>%  filter(Species.code == 'FM'),
  #        aes(x = Log10.mass, y = Log10.P.excretion.rate)) +
  #   geom_point() +
  #   geom_smooth(method = 'lm') +
  #   theme(legend.position = 'none')
  # modelP.FM <- lm(Log10.P.excretion.rate ~ Log10.mass, 
  #                 data = excr %>%  
  #                   filter(Species.code == 'FM'))
  # modelP.FM$coefficients["Log10.mass"]
  
  # ...calculate b coeff of variation for all excretion rates + all species ----
  # to normalize excretion rates
  # create datasets for nutrient and DOM excretion rates separately
  # so that remove specific fish that showed abnormal DOM residuals
  # only for DOM variables
  excr.NPC.var <- excr %>% 
    select(ID, Site.name, 'N.excretion.rate', 'P.excretion.rate', 'C.excretion.rate', Mass)
  # filter out fish species that showed very bad DOM residuals
  excr.DOM.var <- excr %>% 
    select(ID, Site.name, ends_with('excretion.rate'), Mass, 
           -starts_with(c('N.e', 'P.e', 'C.e'))) %>%
    dplyr::filter(
      !is.na(SUVA.excretion.rate),
      # !between(ID, 935, 946),
      # !ID %in% c(1008, 1012, 1013, 115, 901, 904, 907, 911, 
      #            907, 927, 928, 929, 930)
    )
  # check outliers
  coeff_check <- function(x, df) {
    plot <- ggplot(df, aes(x = log10(Mass), y = log10(.data[[x]]))) +
      geom_point() +
      geom_smooth(method = 'lm')
    print(plot)
  }
  for (x in colnames(excr.NPC.var)) {
    if (!x %in% c("ID", "Site.name", "Mass")) {
    coeff_check(x, excr.NPC.var)
    }
  }
  for (x in colnames(excr.DOM.var)) {
    if (!x %in% c("ID", "Site.name", "Mass")) {
      coeff_check(x, excr.DOM.var)
    }
  }
  
  # filter out some outliers
  # excr.DOM.var <- excr.DOM.var %>%
  #   dplyr::filter(log10(C1.excretion.rate) > -4.5)
  # excr.NPC.var <- excr.NPC.var %>% mutate(C.excretion.rate = 
  #                                           case_when(C.excretion.rate ~ 12.5, 
  #                                                   NA, .))
  
  # combine NPC and DOM datasets
  excr.var <- left_join(excr.NPC.var, excr.DOM.var, by = c('ID', 'Mass', 'Site.name'))
  
  # ..loop to generate coefficient of variation for each excretion rate variable ----
  # and use it to mass-normalize each rate in excr.var data frame
  coeff.df <- tibble()
  
  for (col in colnames(excr.var)) {
    # Exclude the ID + Mass columns
    if (!col %in% c("ID", "Site.name", "Mass")) {  
      
      subdf <- excr.var %>% 
        filter(!is.na(.data[[col]]), .data[[col]] > 0)  # Use .data to refer to column
      
      num_observations <- nrow(subdf)
      
      if (num_observations > 11) {
        subdf$log_col <- log10(subdf[[col]])
        subdf$log_Mass <- log10(subdf$Mass)
        
        model <- lm(log_col ~ log_Mass, data = subdf)
        result <- tibble(
          Var = col,
          b.coeff = coef(model)[2],  # Coefficient for log10(Mass)
          n_obs = num_observations
        ) 
        excr.var <- excr.var %>% 
        dplyr::mutate(
          !!paste0(
            "massnorm.", sub("etion.rate", "", col)) := excr.var[[col]]/(Mass^result$b.coeff)
          )
        
      } else {
        result <- tibble(
          Var = col,
          b.coeff = NA,
          n_obs = num_observations
        )
      }
      
      coeff.df <- bind_rows(coeff.df, result)
      cat("Processed column:", col, "\n")
    }
  }
  coeff.df
  excr.var
  
  # calculate mass-normalized ratios
  excr.var <- excr.var %>%
    mutate(
      massnorm.NP.excr = (massnorm.N.excr/massnorm.P.excr)/(14/31),
      massnorm.CN.excr = (massnorm.C.excr/massnorm.N.excr)/(12/14),
      massnorm.CP.excr = (massnorm.C.excr/massnorm.P.excr)/(12/31)
      )
  
  # add newly created columns to excr dataset
  excr <- left_join(excr, excr.var)
  
  ####################### Create datasets for data analysis ###################
  # function to rename some DOM characteristics
  rename_DOM <- function(df) {
    df <- df %>% 
      rename(
        SUVA254 = AmSUVA,
        βα = AmBA
      ) %>% 
      rename_with(~ sub("Am", "", .))
    return(df)
  }
  
  # ..make excr dataset with one entry for each excretion average ----
  # ..N/P excretion species average ----
  excr.sp <- excr %>% 
    group_by(Site.name, Species.code, Trophic.position2, AmDOC, DOC.level) %>% 
    summarise(
      across(
        ends_with('excr'), 
        function(x) median(x, na.rm = TRUE),
        .names = "{.col}.sp"
      ),
      n = n()
      )
  
  # ..pivot dataset for DOM excretion only ----
  excr.DOM <- excr.var %>% 
    select(ID, Site.name, massnorm.SUVA.excr:massnorm.C7.excr) %>% 
    mutate(
      DOC.level = factor(ifelse(Site.name == 'L224', 'low', 
                                ifelse(Site.name == 'L239', 'med',
                                       ifelse(Site.name == 'L222', 'high', NA)))),
      DOC.level = fct_relevel(DOC.level, c('low', 'med', 'high'))
      ) %>% 
    pivot_longer(
      massnorm.SUVA.excr:massnorm.C7.excr,
      names_to = 'type',
      names_pattern = 'massnorm\\.(.*)\\.excr',
      values_to = 'massnorm.excr',
      values_drop_na = TRUE
    ) %>% 
    dplyr::mutate(
      type = if_else(type == 'SUVA', 'SUVA254',
                     if_else(type == 'BA', 'βα', type)),
      type = reorder(type, massnorm.excr, median)
    )
  
  # ..prepare PCA dataset ----
  excr.pca <- excr %>% group_by(Site.name) %>%
    summarise(across(c(
      starts_with('Am'),-AmC6,-starts_with(c(
        'AmA', 'AmP', 'AmS2', 'Ambi', 'AmS3', 'AmHis', 'AmR'
      )),
      Area:Part.P
    ),
    \(x) mean(x, na.rm = TRUE))) 
  
  # transform PARAFAC components into percentages for PCA
  Calc_AmCtot <- function(df) {
    result <- rowSums(df[, c("AmC1", "AmC2", "AmC3", "AmC4", "AmC5", "AmC7")])
    return(result)
  }
  
  for (site in unique(excr.pca$Site.name)) {
    AmCtot <- Calc_AmCtot(excr.pca)
    
    for (col in c("AmC1", "AmC2", "AmC3", "AmC4", "AmC5", "AmC7")) {
      col_name <- paste0(col, "per")
      excr.pca <- excr.pca %>%
        mutate(!!col_name := !!sym(col) / AmCtot * 100)
    }
    excr.pca <- excr.pca %>%
      mutate(C_humicper = AmC1per + AmC2per + AmC3per,
             C_microbialper = AmC5per)
  }
  cat("AmCtot:", AmCtot, "\n")
  
  # ..prepare NMDS dataset ----
  excr.amb <- excr.pca %>% 
    rename(ID = Site.name) %>% 
    filter(ID %in% c('L222', 'L224', 'L239')) %>% 
    dplyr::mutate(Site.name = c('L222', 'L224', 'L239'),
                  Source = c('AmL222', 'AmL224', 'AmL239'),
                  Trophic.position2 = c('AmL222', 'AmL224', 'AmL239')) %>% 
    rename_DOM() %>% 
    select(
      -c(C1per:C_microbialper),
      -c(DOC, TDP, TDN)
    )
  
  excr.nmds <- excr %>% 
    select(c(
      ID, Site.name, Species.code, Trophic.position2,
      ends_with('excr'),
      -ends_with(c('N.excr', 'P.excr', 'C.excr'))
    )
    ) %>%
    mutate(across(where(is.numeric),
                  ~ if_else(. < 0, 0, .)),
           ID = as.character(ID)) %>% 
    rename_with( ~ sub("massnorm.", "", .)) %>%
    rename_with( ~ sub(".excr", "", .)) %>%
    rename(
      Source = Species.code,
      SUVA254 = SUVA,
      βα = BA
    ) %>%
    bind_rows(excr.amb) %>%
    dplyr::filter(!is.na(C4))
   
  # ..summary statistics ----
  excr.ss <- excr %>% 
    select(c('massnorm.N.excr', 'massnorm.P.excr', 'massnorm.NP.excr',
             'Mass')) %>% 
    describe_distribution()
  
  excrtp.ss <- excr %>% group_by(Site.name, Trophic.position2) %>%
    select(c('massnorm.N.excr', 'massnorm.P.excr', 'massnorm.NP.excr',
             'Mass')) %>% 
    describe_distribution()
  