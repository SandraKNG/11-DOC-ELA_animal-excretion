 #### this is a test with an example dataset to leanr how to use staRdom package #####
 
 # load package ----
 library("tidyverse")
 library("staRdom")
 
 # Number of CPU cores to use
 cores <- detectCores(logical = FALSE)
 
 # example PARAFAC model 
 data(pf_models)
 
 # import raw data
 # EEM
 # example
 folder <- system.file("extdata/EEMs/", package = "staRdom") # folder containing example EEMs
 eem_list_ex <- eem_read(folder, recursive = TRUE, import_function = eem_csv) 
 # in case you use your own data, just replace folder by a path. e.g. "C:/folder/another folder" and change import_function according to instrument.
 eem_list <- eem_read("data/EEM_fish_staRdom", 
                      recursive = TRUE, import_function = "cary") 
 
 # ABS
 absorbance_path = system.file("extdata/absorbance", package = "staRdom") 
 absorbance_ex <- absorbance_read(absorbance_path, cores = cores) # load csv or txt tables in folder
 # load example data, set a path without using system.file to use your own data e.g. "C:/folder/another folder"
 absorbance <- absorbance_read("data/ABS_fish_staRdom", 
                                cores = cores)
 absorbance1 <- absorbance_read("data/ABS_fish_staRdom/ELA_JUNE_2022_05_08_2022_staRdom.csv", 
                                cores = cores)
 absorbance2 <- absorbance_read("data/ABS_fish_staRdom/ELA_JUNE_2022_05_08_2022_staRdom2.csv", 
                                cores = cores)
 absorbance3 <- absorbance_read("data/ABS_fish_staRdom/ELA_JUNE_2022_08_08_2022_staRdom.csv", 
                                cores = cores)
 absorbance4 <- absorbance_read("data/ABS_fish_staRdom/ELA_JUNE_2022_09_08_2022_staRdom.csv", 
                                cores = cores)
 absorbance5 <- absorbance_read("data/ABS_fish_staRdom/ELA_JUNE_2022_12_08_2022_staRdom.csv", 
                                cores = cores)
 
 # example of loading a metatable
 metatable <- system.file("extdata/metatable_dreem.csv",package = "staRdom") # path to example data, can be replaced by a path to your own data
 meta <- read.table(metatable, header = TRUE, sep = ",", dec = ".", row.names = 1) # load data
 # to create using your own data
 eem_metatemplate(eem_list, absorbance) %>%
         write.csv(file="data/metatable.csv", row.names = FALSE)
 
 # plot EEM samples
 eem_overview_plot(eem_list_ex, spp = 9, contour = TRUE)
 eem_overview_plot(eem_list, spp = 9, contour = TRUE)
 
 # check for issues in EEM data
 problem <- eem_checkdata(eem_list,absorbance,meta,metacolumns = c("dilution"),error=FALSE)
 
 