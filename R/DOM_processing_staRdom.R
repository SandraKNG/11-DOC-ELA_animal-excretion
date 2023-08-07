  #### running staRdom to analyze DOM characteristics (indices + PARAFAC) 
  #### with ELA data

  # load package ----
  library(tidyverse)
  library(eemR)
  library(staRdom)

  # Number of CPU cores to use
  cores <- detectCores(logical = FALSE)

  # import raw data ----
  # ..EEM
  eem_list <- eem_read("data/DOM-indices_processing/EEM_staRdom",
                       recursive = T, import_function = eem_csv2)
  # plot EEM samples
  # eem_overview_plot(eem_list, spp = 9, contour = TRUE)

  # ..ABS
  absorbance <- absorbance_read("data/DOM-indices_processing/ABS_staRdom", cores = cores)

  # to create one using your own data
  eem_metatemplate(eem_list, absorbance) %>%
    write.csv(file="data/DOM-indices_processing/metatable.csv", row.names = FALSE)
  # to read metatable
  meta <- read.table("data/metatable.csv", header = TRUE, sep = ",",
                     dec = ".", row.names = 1) # load data
  meta <- meta %>% mutate(dilution = 1)

  # check for issues in EEM and ABS data
  problem <- eem_checkdata(eem_list,absorbance,meta,metacolumns = c("dilution"),error=FALSE)

  # data preparation and correction ----
  # 1) remove FD3 from sample EEM names
  eem_list <- eem_name_replace(eem_list,c("\\(FD3\\)"),c(""))

  # 2) Absorbance baseline correction
  absorbance <- abs_blcor(absorbance,wlrange = c(680,700))

  # 3) Spectral correction
  # ..Excitation
  Excor <- data.table::fread("data/DOM-indices_processing/Excorr.csv")
  # ..Emission
  Emcor <- data.table::fread("data/DOM-indices_processing/Emcorr.csv")

  # Adjust EEMs range to cover vector corrections
  eem_list <- eem_range(eem_list,ex = range(Excor[,1]), em = range(Emcor[,1]))
  eem_list <- eem_spectral_cor(eem_list,Excor,Emcor)

  # 4) Blank substraction
  # extending and interpolation data
  eem_list <- eem_extend2largest(eem_list, interpolation = 1, extend = FALSE, cores = cores)

  # blank substraction
  eem_list <- eem_remove_blank(eem_list)
  
  # eem_overview_plot(eem_list, spp=9, contour = TRUE)
  
  # 5) Inner-filter effect (IFE) correction
  #culv = cuvette length (cm)
  eem_list <- eem_ife_correction(eem_list,absorbance, cuvl = 1)
  
  # 6) Raman normalization
  eem_list <- eem_raman_normalisation2(eem_list, blank = "blank")
  # eem_overview_plot(eem_list, spp=9, contour = TRUE)
  
  # Blanks substraction from dataset
  # ..from EEMs
  eem_list <- eem_extract(eem_list, c("nano", "miliq", 
                                      "milliq", "mq", "blank"),ignore_case = TRUE)
  # ..from absorbance
  absorbance <- dplyr::select(absorbance, -matches("nano|
                                                  miliq|]milliq|mq|blank", 
                                                   ignore.case = TRUE))
  # remove and interpolate scattering
  # dispersion removal
  # creation of a vector that indicates if the values we want to remove follow this order:
  # “raman1”, “raman2”, “rayleigh1” and “rayleigh2”
  remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)

  # creation of a vector that indicates the width (nm) of each wavelength to be removed.
  # The order is the same as previously
  remove_scatter_width <- c(15,15,15,15)

  # removal
  eem_list <- eem_rem_scat(eem_list, remove_scatter = remove_scatter,
                           remove_scatter_width = remove_scatter_width)
  # eem_overview_plot(eem_list, spp=9, contour = TRUE)
  
  # Interpolation
  eem_list <- eem_interp(eem_list, cores = cores, type = 1, extend = FALSE)
  
  eem_overview_plot(eem_list, spp=36, contour = TRUE)
  
  # Data smoothing - ONLY for indices, NOT PARAFAC
  # n sets the width of the rolling mean window in nm
  eem4peaks <- eem_smooth(eem_list, n = 4, cores = cores)
  
  # Overview of samples
  summary(eem_list)

  # peak picking and indices ----
  
  # Biological index (BIX), gives infos on autochtonous production/
  # DOM from aquatic microbial sources;
  # 0.6-0.7 = More degraded DOM from terrestrial sources and
  # >1 = freshly produced DOM by aquatic bacteria
  bix <- eem_biological_index(eem4peaks)
  
  # Coble peaks,
  # B = protein-like tyrosine or tryptophan (microbial sources),
  # T = protein-like tryptophan (microbial sources),
  # A = humic substances (vascular plants sources),
  # C = humic substances (vascular plants sources) and
  # M = humic substances (autochtonous production sources)
  coble_peaks <- eem_coble_peaks(eem4peaks)
  
  # Fluorescence index (FI), gives info on DOM sources;
  # 1.7-2 = microbial sources,
  # 1.2-1.5 = soil and terrestrial plants sources
  fi <- eem_fluorescence_index(eem4peaks)
  
  # Humification index (HIX), gives info on the humification degree;
  # 10-16 = terrestrial sources,
  # <4 = autochtonous sources
  hix <- eem_humification_index(eem4peaks, scale = TRUE)
  # creation of the data table
  indices_peaks <- bix %>%
    full_join(coble_peaks, by = "sample") %>%
    full_join(fi, by = "sample") %>%
    full_join(hix, by = "sample")
  
  # data table
  indices_peaks 
  
  # absorbance indices (a254 --> SUVA 254, SR) ----
  # indices calculations
  slope_parms <- abs_parms(absorbance, cuvl = 1, cores = cores) 
  
  # make one dataset for fluorescence and absorbance indices
  DOM_indices <- left_join(indices_peaks, slope_parms, by = 'sample')
  DOM_indices %>% write.csv(file = "data/DOM-indices_processing/DOM_indices.csv",
                            row.names = F)
  
  # PARAFAC ----
  # samples set wavelength ranges (so that same for all samples)
  eem_red2smallest(eem_list)
  
  # find and remove noise in EEMs
  eem_overview_plot(eem_list, contour = T, spp = 36)
  eem_list <- eem_extract(eem_list, c("ELAFISH903", "ELAFISH904", "ELAFISH906",
                                      "ELAFISH907"),ignore_case = TRUE)
  
  # ..1st iteration ----
  # minimum and maximum number of components - between 3 and 7 or more
  dim_min <- 7
  dim_max <- 10
  
  # number of similar models from which the best is picked
  nstart <- 25
  
  # maximum number of iterations within the PARAFAC model
  maxit = 5000
  
  # tolerance of the model in PARAFAC analysis (10^-8 or 10 recommended)
  ctol <- 10^-6
  
  # Models calculation, one for each number of component. 
  # Warnings will pop to increase nstart...will be done in the last steps
  # pf1 <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = FALSE, 
  #                    const = c("uncons", "uncons", "uncons"), maxit = maxit, 
  #                    nstart = nstart, ctol = ctol, cores = cores)
  
  # same model but using non-negative constraints - GOOD ONE
  pf1n <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = FALSE,
                      const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, 
                      nstart = nstart, ctol = ctol, cores = cores)
  
  
  # rescale B and C modes to a maximum fluorescence of 1 for each component
  # pf1 <- lapply(pf1, eempf_rescaleBC, newscale = "Fmax")
  pf1n <- lapply(pf1n, eempf_rescaleBC, newscale = "Fmax")
  
  # This plot is not shown, because the components violate the assumptions 
  # for fluorescence peaks (negative fluorescence). 
  # Please try, if you are interested.
  # eempf_compare(pf1, contour = TRUE)
  
  # good model plot
  eempf_compare(pf1n, contour = TRUE)
  eempf_fits(pf1n[[1]])
  eempf_plot_comps(pf1n)
  
  # check correlation between different components - table
  # here, we choose model 4
  eempf_cortable(pf1n[[4]], normalisation = FALSE)
  eempf_corplot(pf1n[[4]], progress = FALSE, normalisation = FALSE)
  
  # ..2nd iteration ----
  # normalisation = T because of some of the components highly correlated
  pf2 <- eem_parafac(eem_list, comps = seq(dim_min,dim_max), normalise = TRUE, 
                     const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, 
                     nstart = nstart, ctol = ctol, cores = cores)
  
  # rescale B and C modes
  pf2 <- lapply(pf2, eempf_rescaleBC, newscale = "Fmax")
  
  # eempf_compare(pf2, contour = TRUE) # use this to show the same plot as above
  # for now, we are happy with just the components
  eempf_plot_comps(pf2, contour = TRUE, type = 1)
  
  # calculate leverage (outliers) for model 5 with 7 components
  cpl <- eempf_leverage(pf2[[5]])
  
  # plot leverage - displaying outliers
  eempf_leverage_plot(cpl,qlabel=0.1)
  
  # plot leverage, not so nice plot but interactive to select what to exclude
  # saved in exclude, can be used to start over again with 
  # eem_list_ex <- eem_list %>% eem_exclude(exclude) above
  # exclude <- eempf_leverage_ident(cpl,qlabel=0.1)
  
  # samples, excitation and emission wavelengths to exclude, 
  # makes sense after calculation of leverage
  exclude <- list("ex" = c(),
                  "em" = c(),
                  "sample" = c("ELAFISH121","ELAFISH115", "L373JUNE1022", 
                               "ELAFISH926","ELAFISH1003"))
  
  # exclude outliers if neccessary. if so, restart analysis
  eem_list_ex <- eem_exclude(eem_list, exclude)
  
  # ..3rd iteration ----
  # after outliers removal
  pf3 <- eem_parafac(eem_list_ex, comps = seq(dim_min,dim_max), 
                     normalise = TRUE, maxit = maxit, nstart = nstart, 
                     ctol = ctol, cores = cores)
  pf3 <- lapply(pf3, eempf_rescaleBC, newscale = "Fmax")
  
  # plot
  eempf_compare(pf3, contour = TRUE)
  eempf_plot_comps(pf3, contour = TRUE, type = 1)
  
  # check outliers
  eempf_leverage_plot(eempf_leverage(pf3[[3]]),qlabel=0.1)
  
  # examine residuals
  eempf_residuals_plot(pf3[[5]], eem_list, residuals_only = TRUE, 
                       select = c("ELAFISH1008", "L373JUNE2722",
                                  "L470JUNE2022", "ELAFISH118"), spp = 6, 
                       cores = cores, contour = TRUE)
  
  # adding one more to exclusion list (but have to add all the previous ones)
  exclude <- list("ex" = c(),
                  "em" = c(),
                  "sample" = c("ELAFISH121","ELAFISH115", "L373JUNE1022", 
                               "ELAFISH926","ELAFISH1003","L373JUNE2722"))

  # exclude outliers if neccessary. if so, restart analysis
  eem_list_ex <- eem_exclude(eem_list, exclude)
  
  # ..4th iteration ----
  # recalculating the model with increased accuracy 
  # with only model with 6 components
  ctol <- 10^-8 # decrease tolerance in PARAFAC analysis
  nstart = 25 # number of random starts
  maxit = 10000 # increase number of maximum iterations
  
  pf4 <- eem_parafac(eem_list_ex, comps = seq(dim_min,dim_max), normalise = TRUE, 
                     const = c("nonneg", "nonneg", "nonneg"), 
                     maxit = maxit, nstart = nstart, ctol = ctol, 
                     output = "all", cores = cores, strictly_converging = TRUE)
  
  pf4 <- lapply(pf4, eempf_rescaleBC, newscale = "Fmax")
  
  # check convergence behaviour of created models
  eempf_convergence(pf4[[1]])
  
  # just one model, not really a need to compare
  eempf_compare(pf4, contour = TRUE)
  
  # check outliers
  eempf_leverage_plot(eempf_leverage(pf4[[1]])) 
  # [[1]] means the 4th model in the list, 7 component model in that case
  # because there is only 1 type of model here
  
  # check correlation
  eempf_corplot(pf4[[1]], progress = FALSE)
  
  # FINAL PARAFAC MODEL ----
  # plot resulting components and loadings
  eempf_comp_load_plot(pf4[[1]], contour = TRUE)
  # components only
  ggeem(pf4[[1]], contour = TRUE)
  eempf_comps3D(pf4[[1]])
  # loadings only
  eempf_load_plot
  # this function can be used to view the B- and C-modes
  eempf_plot_comps(pf4[1], type = 2) 
  
  # plot components in each sample, residual and whole sample
  eempf_residuals_plot(pf4[[1]], eem_list, select = eem_names(eem_list)[80:84], 
                       cores = cores, contour = TRUE)
  
  # ..model validation ----
  # calculate split_half analysis to show stability of model
  # data is recombined in 6 different ways and results from each sub-sample
  # should be similar
  sh <- splithalf(eem_list_ex, 7, normalise = TRUE, rand = FALSE, cores = cores,
                  nstart = nstart, strictly_converging = TRUE, maxit = maxit, 
                  ctol = ctol)
  # plot results
  # model is stable if graphs of all components look similar
  splithalf_plot(sh_r)
  
  # check Tucker's Congruency Coefficients: 1 would be perfect similarity
  tcc_sh_table <- splithalf_tcc(sh_r)
  tcc_sh_table
  
  # only if some splits look different, change rand = F to T
  sh_r <- splithalf(eem_list_ex, 5, normalise = TRUE, rand = TRUE, cores = cores, 
                    nstart = nstart, maxit = maxit, ctol = ctol)
  
  # importance of components
  varimp <- eempf_varimp(pf4[[1]], eem_list_ex, cores = cores)
  varimp
  
  # NAMING MODEL + PARAFAC ----
  # get current model names (none set so far!)
  names(pf3)
  # set new model names, number of models must be equal to number of names
  names(pf4) <- c("7 components")
  names(pf4)
  
  # get current component names
  eempf_comp_names(pf4)
  # set new model components names, 
  # renaming the different components name in the final model
  eempf_comp_names(pf4) <- c("E1","E2","E3","E4", "E5", "E6", "E7")
  
  # renaming the components name in several models at a time (ex: pf3)
  eempf_comp_names(pf3) <- list(c("E1","E2","E3","E4", "E5", "E6", "E7"), # model 1
                                c("humic","T2","whatever","peak"),
                                c("rose","peter","frank","dwight","susan"),
                                c("A4","B4","C4","D4","E4","F4"),
                                c("A5","B5","C5","D5","E5","F5","G5"))
  
  # Final model plot with new names
  pf4[[1]] %>%
    ggeem(contour = TRUE)
  
  # components can be reordered
  eempf_reorder
  
  # EXPORTING + INTERPRETING MODEL ----
  # comparing your results with the openfluor.org database
  eempf_openfluor(pf4[[1]], 
                  file = "data/DOM-indices_processing/PARAFAC_staRdom/7components_model_openfluor.txt")
  
  # exporting model matrices with samples to csv file
  eempf_export(pf4[[1]], 
               export = "data/DOM-indices_processing/PARAFAC_staRdom/7components_model_openfluor.csv")
  
  # export final samples results with metdadata table
  eem_metatemplate