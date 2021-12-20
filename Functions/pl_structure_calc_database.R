
# 1. Load necessary functions ---------------------------------------------

source(file = here("Functions","Functions_for_pl_structure_calc.R"))

# Fri Jan 11 16:31:07 2019 ------------------------------

# 2. pl_structure_calc function definition --------------------------------


###################################################################################
#                                                                                 #        
#"#"structure calc funciton, extraction MS2 information for the s --------        #
#                                                                                 #
###################################################################################


pl_structure_calc <- function(ms2_spectraidentifier, range_carb, range_double,samplename){
  
  # source(file = here("scripts","Functions_for_pl_structure_calc.R"))
  # 1. select the mass correction values -------------
  
  ## m/z calibration according to the fragments of the internal standard
  rm(cal_spectra)
  cal_spectra <- ms2_spectraidentifier %>% filter(Name %in% c("PC:36:2:0:3","PE:36:2:0:0","PC:28:0:0:0","PE:28:0:0:0","PEO:36:2:0:0","PEO:36:3:0:0"))
  cal_spectra$seqNum
  # Collect all internal standard spectra
  cal_ref <- cal_obs <- c() # Initiate calibration vectors
  
  if(nrow(cal_spectra)>=2){ # Checks if an internal standard calibration spectrum was recorded
    cal_ref <- c(227.2011,255.2324,277.2168,279.2324, 281.2481,283.2637,307.2637,327.2324,462.2985,464.3141) # Exact masses of three reference MS/MS fragment peaks from the internal standard 
    cal_spectra$MS2 <- peaks(object = mz,scans= cal_spectra$seqNum)
    cal_obs <- cal_spectra %>% mutate(obs=future_map(MS2,select_nearest_peak)) # Extracts actually measured m/z of reference peaks
    
  }else{
    print("ERRRRRRORRRR")
    # cal_ref <- c(591.4026, 647.4288, 727.3951) # Select three reference MS/MS fragment peaks from internal standard 
    # cal_obs <- c(591.4600, 647.55, 727.50) # For the case no internal standard was used, add alternative calibration data here
  }
  cal_obs <- bind_rows(cal_obs$obs)
  fit <- tibble(y=cal_obs$cal_p,x=cal_obs$cal_o)
  cal_fit <<- lm(y ~ poly(x, 1, raw=TRUE), data = fit) # Generates calibration model based on measured and theoretical fragment m/z
  #summary(cal_fit)
  # plot(x = fit$x,y = fit$y)
  
  # 2. MS2 daten als List column angehängt. -----------
  ms2_spectraidentifier$MS2_data <- peaks(object = mz,scans= ms2_spectraidentifier$seqNum)
  
  # 3. Anwenden der masscorrections funktion, originale beleiben erhalten -------
  ms2_spectraidentifier <-   ms2_spectraidentifier %>% mutate(MS2_data_masscorr=future_map(MS2_data,correct_massdrift))
  rm(cal_fit)
  # 4. erstellen aller möglichen FS seitenketten Kombinationen und dem FA_mass_dataframe --------
  n.ms2_spectraidentifier <- ms2_spectraidentifier 
  #data preparation for FA.comb function
  n.ms2_spectraidentifier.2 <- n.ms2_spectraidentifier %>% mutate(realC=str_extract(string = Name,pattern = "(?<=:)[0-9]{2}"),realDB=str_extract(string = Name,pattern = "(?<=:[0-9]{2}:)[0-9]{1,2}"))
  n.ms2_spectraidentifier.3 <- n.ms2_spectraidentifier.2 %>% mutate(Name1=Name) %>% group_by(Name1) %>% nest()
  # Thu Nov 15 12:02:49 2018 ------------------------------
  sql.export <- c()
  load(file = here("database",paste0("Run",Sys.Date(),".Rdata")))
  sql.export_add <- n.ms2_spectraidentifier %>% select(-MS2_data)
  sql.export <- bind_rows(sql.export,sql.export_add)
  save(x = sql.export ,file = here("database",paste0("Run",Sys.Date(),".Rdata")))
}
