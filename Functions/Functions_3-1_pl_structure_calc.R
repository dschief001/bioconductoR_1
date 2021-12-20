# Fri Jan 11 16:31:07 2019 ------------------------------
# Wed Nov 27 17:47:33 2019 ------------------------------

# 2. pl_structure_calc function definition --------------------------------

###################################################################################
#                                                                                 #        
#"#"structure calc funciton, extraction MS2 information for the s --------        #
#                                                                                 #
###################################################################################

pl_structure_calc <- function(ms2_spectraidentifier, range_carb, range_double,samplename,Date_today){
  
  source(file = here("Functions","Functions_3-2_for_pl_structure_calc.R"))
  # 1. select the mass correction values -------------

  ## m/z calibration according to the fragments of the internal standard
  rm(cal_spectra)
  # cal_spectra <- ms2_spectraidentifier %>% filter(Name %in% c("PC:36:2:0:3","PE:36:2:0:0","PC:28:0:0:0","PE:28:0:0:0","PEO:36:2:0:0","PEO:36:3:0:0"))
  # cal_spectra <- ms2_spectraidentifier %>% filter(Name %in% c("PC 36:2","PE 36:2","PC 28:0","PE 28:0","PEO 36:2","PEO 36:3"))
  # cal_spectra <- ms2_spectraidentifier %>% filter(Name %in% c("PC:36:2","PE:36:2","PC:28:0","PE 28:0","PEO:36:2","PEO:36:3"))
  cal_spectra <- ms2_spectraidentifier %>% filter(Name %in% c("PC:36:2","PE:36:2","PC:28:0","PE:28:0","PC:36:4","PE:36:2","PE:38:6","PC:38:6","PC:34:1","PE:34:1","PEO:36:2","PEP:36:2"))
  # cal_spectra$seqNum
  # Collect all internal standard spectra
  cal_ref <- cal_obs <- c() # Initiate calibration vectors

  if(nrow(cal_spectra)>=2){ # Checks if an internal standard calibration spectrum was recorded
    cal_ref <- c(227.2011,255.2324,277.2168,279.2324, 281.2481,283.2637,307.2637,327.2324,462.2985,464.3141) # Exact masses of three reference MS/MS fragment peaks from the internal standard
    cal_spectra$MS2 <- peaks(object = mz,scans= cal_spectra$seqNum)
    cal_obs <- cal_spectra %>% mutate(obs= map(MS2,select_nearest_peak)) # Extracts actually measured m/z of reference peaks
    cal_obs <- bind_rows(cal_obs$obs)
    fit <- tibble(y=cal_obs$cal_p,x=cal_obs$cal_o)
  }else{
    print("calibration ERRRRRRORRRR")
    cal_ref <- c(227.2011,255.2324,277.2168,279.2324, 281.2481,283.2637,307.2637,327.2324,462.2985,464.3141) # Exact masses of three reference MS/MS fragment peaks from the internal standard
    cal_obs <- c(227.2011,255.2324,277.2168,279.2324, 281.2481,283.2637,307.2637,327.2324,462.2985,464.3141) # Exact masses of three reference MS/MS fragment peaks from the internal standard
    # cal_ref <- c(591.4026, 647.4288, 727.3951) # Select three reference MS/MS fragment peaks from internal standard
    # cal_obs <- c(591.4600, 647.55, 727.50) # For the case no internal standard was used, add alternative calibration data here
    fit <- tibble(y=cal_ref,x=cal_obs)
  }
  cal_fit <<- lm(y ~ poly(x, 1, raw=TRUE), data = fit) # Generates calibration model based on measured and theoretical fragment m/z
  #summary(cal_fit)
  # plot(x = fit$x,y = fit$y)

  # 2. MS2 daten als List column angehängt. -----------
  ms2_spectraidentifier$MS2_data <- peaks(object = mz,scans= ms2_spectraidentifier$seqNum)
  # 3. Anwenden der masscorrections funktion, originale beleiben erhalten -------
  ms2_spectraidentifier <-   ms2_spectraidentifier %>% mutate(MS2_data_masscorr= map(MS2_data,correct_massdrift))
  rm(cal_fit)
  
  print("===                                                       8%")
  # 4. erstellen aller möglichen FS seitenketten Kombinationen und dem FA_mass_dataframe --------
  #data preparation for FA.comb function
  
  # implement the vinyl DB here:
  n.ms2_spectraidentifier.2 <- ms2_spectraidentifier %>% mutate(
    realC=str_extract(string = Name,pattern = "(?<=:)[0-9]{2}"),
    realDB=if_else(condition = str_sub(string = ms2_spectraidentifier$Name,start = 3,end = 3)=="P",
                   true = as.numeric(str_extract(string = Name,pattern = "(?<=:[0-9]{2}:)[0-9]{1,2}"))+1,
                   false = as.numeric(str_extract(string = Name,pattern = "(?<=:[0-9]{2}:)[0-9]{1,2}"))
    ))
  n.ms2_spectraidentifier.2 <- n.ms2_spectraidentifier.2 %>% mutate(Name1=Name) %>% group_by(Name1) %>% nest()
  # Thu Nov 15 12:02:49 2018 ------------------------------
  # need to change here the possible FA combination for PE[e/p]
  
  # TIME consuming step 1 ------------------------------------------------
  n.ms2_spectraidentifier.2 <- n.ms2_spectraidentifier.2 %>% mutate(FA_possible= map(.x = data,.f = generate_possible_FA_layout))
  n.ms2_spectraidentifier.2 <- n.ms2_spectraidentifier.2 %>% unnest(cols = data)
  # Initialize theoretically possible fragment structure space using the function fragments_theo()
  
  # TIME consuming step 2 ------------------------------------------------
  print("=======                                                  20%")
  
  n.ms2_spectraidentifier.2 <- n.ms2_spectraidentifier.2 %>% mutate(FA_Frag_mass= list(fragments_theo(carb = range_carb, double = range_double)),Ether_restmass=list(fragments_theo_ether(carb = range_carb, double = range_double)))
  
  n.ms2_spectraidentifier.2 <- n.ms2_spectraidentifier.2 %>% group_by(Name) %>% nest()
  # ! check for intensity used! 5. Filter only the possible FA combinations     -----------------------------
  n.ms2_spectraidentifier.7 <- n.ms2_spectraidentifier.2 %>% mutate(FA_Frag_massN= map(data, actual_filterfunction),M_Frag_massN= map(data,actual_filterfunction_ether))
  rm("n.ms2_spectraidentifier.2")
  #=> should store in QC log this new file
  n.ms2_spectraidentifier.7 <- n.ms2_spectraidentifier.7 %>% unnest(cols=data)
  
  n.ms2_spectraidentifier.7 <- n.ms2_spectraidentifier.7 %>% select(Name,seqNum,totIonCurrent,precursorMZ,FA_possible,MS2_data_masscorr,FA_Frag_massN,M_Frag_massN,FA_Frag_mass,Ether_restmass)
  
  print("============                                             30%")
  
  ##### Overlaps, double assigned peaks are now pooled
  # !!! double assigned peaks are found in n.ms2_spectraidentifier.8 --------
  n.ms2_spectraidentifier.8use <- n.ms2_spectraidentifier.7 %>% group_by(Name1_=Name,seqNum) %>% nest()
  # 6. Actual FA to FA_mass_table (not combined) --------
  n.ms2_spectraidentifier.9 <- n.ms2_spectraidentifier.8use %>% mutate(Found_FAs= map(data,ms2_annotation))
  
  # TIME Feedback plots -----------------------------------------------------
  Feedback_plots(n.ms2_spectraidentifier.9)
  Contribution <<- get_contribution_data(n.ms2_spectraidentifier.9)
  # TIME Contribution plots -----------------------------------------------------
  Contribution_plot(Contribution)

  n.ms2_spectraidentifier.10 <- n.ms2_spectraidentifier.9 %>% unnest(cols=data) # ,.preserve = Found_FAs)
  n.ms2_spectraidentifier.11 <- n.ms2_spectraidentifier.10 %>% ungroup() %>% select(Name,seqNum,totIonCurrent,FA_possible, Found_FAs,FA_Frag_massN,M_Frag_massN)
  # umsortieren der MS2_FA_possible und duplikate entferenn -----------------

  print("===================                                      40%")
  # 7. Cleanup der DAten ----------------------------------------------------
  # !!! careful sortfunction on text does not distinguish between 1  --------
  n.ms2_spectraidentifier.11$FA_possible <- lapply(n.ms2_spectraidentifier.11$FA_possible, cleanup_FA_possible)
  # 8.sort by name ---------------------------------------------------------
  n.ms2_spectraidentifier.11 <- n.ms2_spectraidentifier.11 %>% arrange(Name)

  # 9. generate emty table for data finding together -----------------------
  n.ms2_spectraidentifier.11$FA_possible2 <- lapply(n.ms2_spectraidentifier.11$FA_possible,FUN = function(x) {
    matrix(data = 0, nrow = dim(x)[1], ncol = dim(x)[2])})
  # 10. FA_zugeordnet in die FA-FA_combinations schieben und normali --------
  # zuerst diagnostic data wieder rausholen:
  n.ms2_spectraidentifier.12 <- n.ms2_spectraidentifier.11 %>% mutate(Found_FAs= map(Found_FAs,.f = remove_diagdata))
  n.ms2_spectraidentifier.13 <- n.ms2_spectraidentifier.12
  
  n.ms2_spectraidentifier.13$LAST_FA_profiles_intens <- pmap(.l = list(
    Name=n.ms2_spectraidentifier.12$Name,
    a=n.ms2_spectraidentifier.12$FA_possible,
    b=n.ms2_spectraidentifier.12$FA_possible2,
    c= n.ms2_spectraidentifier.12$Found_FAs,
    e=n.ms2_spectraidentifier.12$totIonCurrent),.f = sort_FAS_into_emtytable_multiplyPrecIntens)

  n.ms2_spectraidentifier.14 <- n.ms2_spectraidentifier.13

  # create new column with plots included.
  print("=============================                            60%")
  n.ms2_spectraidentifier.14$plotlist <-  pmap(list(data=n.ms2_spectraidentifier.13$LAST_FA_profiles_intens,
                                                          seqNum = n.ms2_spectraidentifier.13$seqNum,
                                                          Name = n.ms2_spectraidentifier.13$Name),.f =plot_FA_results )
  # TIME Plot creation ------------------------------------------------------

  # separate plots by type and plot
      # tictoc::tic()
      # plotting_function(data = n.ms2_spectraidentifier.14)
      # tictoc::toc()

  # 11. Abspecheichern des ZwischenREsultats in eine Datei pro Spekt --------
  ms2_saved_datafile <- n.ms2_spectraidentifier.14

  if (dir.exists(path = here("output",paste("MS2_dataextraction_",Date_today),samplename))) {
      save(ms2_saved_datafile,file = here("output",paste("MS2_dataextraction_",Date_today),samplename,paste0(samplename,"_RAW_MS2_dataextraction.RData")))
    }else{
      dir.create(path = here("output",paste("MS2_dataextraction_",Date_today)))
      dir.create(path = here("output",paste("MS2_dataextraction_",Date_today),samplename))
      # dir.create(path = here("output",paste("MS2_dataextraction_",Date_today)))
      save(ms2_saved_datafile,file = here("output",paste("MS2_dataextraction_",Date_today),samplename,paste0(samplename,"_RAW_MS2_dataextraction.RData")))
    }
  print("==========================================               80%")
  # 12. Gewichtet über die MS2 pro Peak mitteln -----------------------------
  n.ms2_F.1 <- n.ms2_spectraidentifier.13 %>% select(Name,seqNum,LAST_FA_profiles_intens)
  n.ms2_F.1 <- n.ms2_F.1 %>% mutate(LAST_FA_profiles_intens2= map(.x = LAST_FA_profiles_intens,.f = filter_average_MS2_data))
  n.ms2_F.2 <- n.ms2_F.1 %>% select(-c(LAST_FA_profiles_intens,seqNum))
  n.ms2_F.2 <- n.ms2_F.2 %>% group_by(Name) %>% nest()
  n.ms2_F.3 <- n.ms2_F.2
  n.ms2_F.3$Final <-  map2(.x = n.ms2_F.2$Name,.y = n.ms2_F.2$data,.f = averaging_over_spectra)
  n.ms2_F.4 <- n.ms2_F.3 %>% select(-data) %>% unnest(cols = Final)
  n.ms2_F.4$position[is.na(n.ms2_F.4$position)] <- "both"
  plotting_function2(n.ms2_F.4)
  
  # Plasmalogen retreatment:
  n.ms2_F.4 <- bind_rows(
    n.ms2_F.4 %>% filter(!(str_detect(string = Name,pattern = "(?<=[A-Z]{2})[pP]{1}") & position=="sn1_FA")),
    n.ms2_F.4 %>%
      filter(
        str_detect(string = Name, pattern = "(?<=[A-Z]{2})[pP]{1}"),
        position == "sn1_FA",
        FA_norm != 0,
      ) %>% 
      filter(str_remove(string = FA, pattern = ".*:") >= 1) %>% mutate(FA =
                                                                         paste0(
                                                                           str_remove(string = FA, pattern = ":.*"),
                                                                           ":",
                                                                           as.numeric(str_remove(string = FA, pattern = ".*:")) - 1
                                                                         ))
  ) %>% group_by(Name) %>% mutate(FA_norm = FA_norm/sum(FA_norm))
    
    
  # 14. Ausgabe des short REsultats. ----------------------------------------
  gc()
  print("========================================================= 100%")
  return(n.ms2_F.4)
}





# Function for MS2 readout combination and distribution of the data -------

# for (i in seq_along(MS2.0$data)) {
#   x <- MS2.0$data[[i]]
#   clean_summariseMS2(x)
#   print(i)
# }

clean_summariseMS2 <- function(x) {
  if (all(is.na(x$FA_norm))) {
    tibble(
      FA = "unknown",
      position = "both",
      FA_Area = unique(x$Area),
      FA_pmol = unique(x$pmol)
    )
  }
  else if (any(is.nan(x$FA_norm))) {
    stop("Error NaN and values mixed in one lipid!")
  }
  else {
    if (all(unique(x$position) %in% c("sn1_FA", "sn2_FA"))) {
      x %>% filter(FA_norm != 0) %>% mutate(
        FA_norm = FA_norm / sum(FA_norm),
        FA_Area = Area * FA_norm,
        FA_pmol = FA_norm * pmol
      ) %>% select(-Area, -pmol, -FA_norm)
    }
    else if (unique(x$position) == "both") {
      x %>% filter(FA_norm != 0) %>%
        mutate(FA_norm = FA_norm / sum(FA_norm),FA_Area = Area * FA_norm, FA_pmol =FA_norm * pmol) %>%
        select(-Area, -pmol, -FA_norm)
    }
  }
}