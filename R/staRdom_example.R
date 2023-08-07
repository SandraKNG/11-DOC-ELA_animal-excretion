 #### this is a test with an example dataset + my own ELA dataset
 ### to learn how to use staRdom package #####
 
 # load package ----
 library(tidyverse)
 library(eemR)
 library(staRdom)
 
 # Number of CPU cores to use
 cores <- detectCores(logical = FALSE)
 
 # example PARAFAC model 
 data(pf_models)
 
 # import raw data ----
 # ..EEM ----
 # example
 folder <- system.file("extdata/EEMs/", package = "staRdom") # folder containing example EEMs
 eem_list_ex <- eem_read(folder, recursive = TRUE, import_function = eem_csv) 
 # plot EEM samples
 eem_overview_plot(eem_list_ex, spp = 9, contour = TRUE)
 
 # ..ABS ----
 # example
 absorbance_path = system.file("extdata/absorbance", package = "staRdom") 
 absorbance_ex <- absorbance_read(absorbance_path, cores = cores) # load csv or txt tables in folder
 
 # loading a metatable
 # example of loading a metatable
 metatable <- system.file("extdata/metatable_dreem.csv",package = "staRdom") # path to example data, can be replaced by a path to your own data
 meta_ex <- read.table(metatable, header = TRUE, sep = ",", dec = ".", row.names = 1) # load data
 
 # check for issues in EEM and ABS data
 problem_ex <- eem_checkdata(eem_list_ex,absorbance_ex,meta_ex,metacolumns = c("dilution"),error=FALSE)
 
 # data preparation and correction ----
 # remove FD3 from EEM names
 # example
 eem_list_ex <- eem_name_replace(eem_list_ex,c("\\(FD3\\)"),c("")) 
 
 # Absorbance baseline correction
 # example
 absorbance_ex <- abs_blcor(absorbance_ex,wlrange = c(680,700))

 # Spectral correction ----
 # example
 # ..Excitation
 excorfile <- system.file("extdata/CorrectionFiles/xc06se06n.csv",package="staRdom")# path
 Excor <- data.table::fread(excorfile)#data
 # ..Emission
 emcorfile <- system.file("extdata/CorrectionFiles/mcorrs_4nm.csv",package="staRdom")#path
 Emcor <- data.table::fread(emcorfile)#data

 
 # Adjust EEMs range to cover vector corrections
 # example
 eem_list_ex <- eem_range(eem_list_ex,ex = range(Excor[,1]), em = range(Emcor[,1]))#create spectral range
 eem_list_ex <- eem_spectral_cor(eem_list_ex,Excor,Emcor)#correction

 
 # Blank substraction ----
 # extending and interpolation data
 eem_list_ex <- eem_extend2largest(eem_list_ex, interpolation = 1,
                                   extend = FALSE, cores = cores)

 # blank substraction
 eem_list_ex <- eem_remove_blank(eem_list_ex)
 eem_overview_plot(eem_list_ex, spp=9, contour = TRUE)
 
 # inner-filter effect (IFE) correction
 eem_list_ex <- eem_ife_correction(eem_list_ex,absorbance_ex, cuvl = 5)#culv = cuvette length (cm)

 # Raman normalization
 eem_list_ex <- eem_raman_normalisation2(eem_list_ex, blank = "blank")
 #Raman correction, blank = correction method (here with the blanks)

 # Blanks substraction
 #from EEMs
 eem_list_ex <- eem_extract(eem_list_ex, c("nano", "miliq", "milliq",
                                           "mq", "blank"),ignore_case = TRUE)
 #from absorbances
 absorbance_ex <- dplyr::select(absorbance_ex, -matches("nano|miliq|milliq|mq|blank",
                                                        ignore.case = TRUE))
 # remove and interpolate scattering
 #dispersion removal

 #creation of a vector that indicates if the values we want to remove follow this order:
 #“raman1”, “raman2”, “rayleigh1” and “rayleigh2”
 remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)

 #creation of a vector that indicates the width (nm) of each wavelength to be removed.
 #The order is the same as previously
 remove_scatter_width <- c(15,15,15,15)

 #removal
 eem_list_ex <- eem_rem_scat(eem_list_ex, remove_scatter = remove_scatter,
                             remove_scatter_width = remove_scatter_width)

 eem_overview_plot(eem_list_ex, spp=9, contour = TRUE)

 #Interpolation

 eem_list_ex <- eem_interp(eem_list_ex, cores = cores, type = 1, extend = FALSE)
 eem_overview_plot(eem_list_ex, spp=9, contour = TRUE)

 # dilution correction
 dil_data <- meta_ex["dilution"]# creation of the dataset containing the dilution factors
 eem_list_ex <- eem_dilution(eem_list_ex,dil_data)#the correction
 eem_overview_plot(eem_list_ex, dil_data)

 # overview of data before indices
 summary(eem_list_ex)
 
 # PARAFAC ----
 # example
 # Name creation for temporary file
 dreem_raw <- tempfile()
 
 download.file("https://gitlab.com/dreem/drEEM/-/raw/master/tutorials_demos/datasets/PortSurveyData_corrected.mat?inline=false",dreem_raw)
 
 dreem_data <- dreem_raw %>%
         R.matlab::readMat()
 
 # erase temporary path to the data since we do not need it anymore
 unlink(dreem_raw)
 
 # creation of a eemlist object
 eem_list <- lapply(dreem_data$filelist.eem, function(file){
         n <- which(dreem_data$filelist.eem == file)
         file <- file %>%
                 gsub("^\\s+|\\s+$", "", .) %>% # erase space in the file names
                 sub(pattern = "(.*)\\..*$", replacement = "\\1", .) # erase extensions from file names
         eem <- list(file = paste0("drEEM/dataset/",file),sample = file,x = dreem_data$XcRU[n,,] %>%
                             as.matrix(),ex = dreem_data$Ex %>% as.vector(), em = dreem_data$Em.in %>%
                             as.vector(), location = "drEEM/dataset/")
         class(eem) <- "eem"
         attr(eem, "is_blank_corrected") <- TRUE
         attr(eem, "is_scatter_corrected") <- FALSE
         attr(eem, "is_ife_corrected") <- TRUE
         attr(eem, "is_raman_normalized") <- TRUE
         attr(eem, "manufacturer") <- "unknown"
         eem
 }) %>%
         `class<-`("eemlist")
 
 # add a prefix "d" to the file names, because R doesn't like when they start with a number
 eem_names(eem_list) <- paste0("d",eem_names(eem_list))
 
 #For this tutorial with drEEM dataset, we need to remove samples containing "bl" and  "0A" in their name
 ol <- function(x){x==("bl") | x == "0A"}
 extract <- dreem_data$sites %>% unlist() %>% ol() %>% which()
 eem_list <- eem_list %>% eem_extract(extract)
 
 # my own data
 # load corrected EEMs
 eem_data <- R.matlab::readMat('data/DOM-indices_processing/PARAFAC/
                               renamed_corrected_EEMs/corrected_EEMS.mat')
 
 eem_list <- lapply(eem_data$filelist.eem, function(file){
         n <- which(eem_data$filelist.eem == file)
         file <- file %>%
                 gsub("^\\s+|\\s+$", "", .) %>% # erase space in the file names
                 sub(pattern = "(.*)\\..*$", replacement = "\\1", .) # erase extensions from file names
         eem <- list(file = paste0("data/DOM-indices_processing/PARAFAC/",file),sample = file,x = eem_data$XcRU[n,,] %>%
                             as.matrix(),ex = eem_data$Ex %>% as.vector(), em = eem_data$Em.in %>%
                             as.vector(), location = "data/DOM-indices_processing/PARAFAC/")
         class(eem) <- "eem"
         attr(eem, "is_blank_corrected") <- TRUE
         attr(eem, "is_scatter_corrected") <- FALSE
         attr(eem, "is_ife_corrected") <- TRUE
         attr(eem, "is_raman_normalized") <- TRUE
         attr(eem, "manufacturer") <- "unknown"
         eem
 }) %>%
         `class<-`("eemlist")
 
 eem_list <- eem_read("data/DOM-indices_processing/PARAFAC/corrected_EEMs", 
                      recursive = F, import_function = 'cary') 
  
 
 