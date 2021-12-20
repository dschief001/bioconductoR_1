source(file = here("Functions","Functions_for_pl_structure_calc.R"))
# p) ms2_annotation_with guessing part only -------------------------------

ms2_annotation_IDonly <- function(data){
   cont_function <- function(identificator,Diff){
      temp <- (identificator %>% mutate(Test=abs(V1-(data$precursorMZ - Diff))) %>%
         filter( mass_tol_identification >= Test) %>% filter(Test==min(Test)) %>%
         select(-Test))$V2/sum(identificator$V2)
      if(length(temp)>0){temp}else{0}
    }
    count_function <- function(identificator,Diff){
      temp <- (identificator %>% mutate(Test=abs(V1-(data$precursorMZ - Diff))) %>%
         filter( mass_tol_identification >= Test) %>% filter(Test==min(Test)) %>%
         select(-Test))$V2
      if(length(temp)>0){temp}else{0}
    }
  if (nrow(data$MS2_data[1][[1]])!=0) {
  measured_masses_gethere <- as_tibble(data$MS2_data[1][[1]])
  # negativ beispiel
  # Peparatatory guess Lipidclass from data
  identificator <- measured_masses_gethere %>% filter(!(round(data$precursorMZ,digits= 0) == round(measured_masses_gethere$V1,digits = 0)) & V1 > lower_threshold_Identification)
  # min(abs(identificator$V1-(data$precursorMZ - 60)))== abs(identificator$V1-(data$precursorMZ - 60))
   
    
    Find_PC_PS <- tibble(
      PC_cont= cont_function(identificator = identificator,Diff = 60),
      PS_cont= cont_function(identificator = identificator,Diff = 87),
      PC_count= count_function(identificator = identificator,Diff = 60),
      PS_count=count_function(identificator = identificator,Diff = 87))
  }else{
    
    Find_PC_PS <-tibble(PC_cont=0,PS_cont=0,
                                   PC_count=0,
                                 PS_count=0)
  }
  return(Find_PC_PS)
}

# y <- list()
# i=1
# for (i in seq_along(n.ms2_spectraidentifier.8use$data)) {
#   print(i)
#   y[[i]] <- local(
#     ms2_annotation_IDonly(data = n.ms2_spectraidentifier.8use$data[[i]])
#   )
#   
# }


# 2. pl_structure_calc function definition --------------------------------


###################################################################################
#                                                                                 #        
#"#"structure calc funciton, extraction MS2 information for the s --------        #
#                                                                                 #
###################################################################################


pl_structure_calc <- function(ms2_spectraidentifier, range_carb, range_double,samplename){
  
  source(file = here("Functions","Functions_for_pl_structure_calc.R"))
  # 1. select the mass correction values -------------
  
  ## m/z calibration according to the fragments of the internal standard
  rm(cal_spectra)
  # cal_spectra <- ms2_spectraidentifier %>% filter(Name %in% c("PC:36:2:0:3","PE:36:2:0:0","PC:28:0:0:0","PE:28:0:0:0","PEO:36:2:0:0","PEO:36:3:0:0"))
  cal_spectra <- ms2_spectraidentifier %>% filter(round(precursorMZ,digits = 1) %in% c(678.5,830.6,742.6,634.5,722.5))
  # %>% group_by(round(precursorMZ,digits = 1)) %>% filter(precursorIntensity==max(precursorIntensity))
  # this line offers the option to filter for only the most intense spectrum found at the given mass.
  head(ms2_spectraidentifier$precursorMZ)
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

  # 4. Save of MS2data into database. --------
  # Thu Nov 15 12:02:49 2018 ------------------------------
  sql.export <- c()
  
  tryCatch({
    load(file = here("database",paste0("Run-total-guessing",Sys.Date(),".Rdata")))
  }, error = function(e) {
    sql.export <- tibble()
  })
  sql.export_add <- ms2_spectraidentifier %>% select(-MS2_data)
  sql.export <- bind_rows(sql.export,sql.export_add)
  save(x = sql.export ,file = here("database",paste0("Run-total-guessing",Sys.Date(),".Rdata")))
  
  # 2. MS2 daten als List column angehängt. -----------
  ms2_spectraidentifier$MS2_data <- peaks(object = mz,scans= ms2_spectraidentifier$seqNum)
  tryCatch({
    load(file = here("database",paste0("Run-allMS2",Sys.Date(),".Rdata")))
  }, error = function(e) {
    sql.export <- tibble()
  })
  sql.export_add <- ms2_spectraidentifier
  sql.export <- bind_rows(sql.export,sql.export_add)
  save(x = sql.export ,file = here("database",paste0("Run-allMS2",Sys.Date(),".Rdata")))

  n.ms2_spectraidentifier.7 <- ms2_spectraidentifier 
  # %>% select(Sample,seqNum,precursorMZ,MS2_data,totIonCurrent,peaksCount,retentionTime,basePeakMZ,basePeakIntensity)
  ##### Overlaps, double assigned peaks are now pooled
  # !!! double assigned peaks are found in n.ms2_spectraidentifier.8 --------
  n.ms2_spectraidentifier.8use <- n.ms2_spectraidentifier.7 %>% group_by(seqNum1=seqNum) %>% nest()
  # 6. Actual FA to FA_mass_table (not combined) --------
  n.ms2_spectraidentifier.9 <- n.ms2_spectraidentifier.8use %>% mutate(Found_Peaks = future_map(data,ms2_annotation_IDonly))
  n.ms2_spectraidentifier.10 <- n.ms2_spectraidentifier.9 %>% unnest() %>% select(retentionTime,precursorMZ,lowMZ, highMZ ,PC_cont,PS_cont , PC_count , PS_count )
  n.ms2_spectraidentifier.10 <- n.ms2_spectraidentifier.10 %>% mutate(retentionTime=retentionTime/60)

  n.ms2_spectraidentifier.11 <- n.ms2_spectraidentifier.10 %>% filter(PC_cont>=0.5) %>% select(Mass=precursorMZ,Time=retentionTime) %>% mutate(Mass=round(Mass,digits = 2),Time=round(x = Time,digits = 2),Name=c(paste0("PC_",seq_along((n.ms2_spectraidentifier.10 %>% filter(PC_cont>=0.5))$retentionTime))),Annotations="MS2derived")
  # temp <- n.ms2_spectraidentifier.11 %>% mutate(PS=as.numeric(str_remove(string = n.ms2_spectraidentifier.11$Mass,pattern = "\\..*"))%%2)
  if (!dir.exists(paths = here("Output",paste0("MS2_dataextraction_",Sys.Date(),"_Peaklists")))){
    dir.create(path = here("Output",paste0("MS2_dataextraction_",Sys.Date(),"_Peaklists")))  }
  write.csv(n.ms2_spectraidentifier.11,file = here("Output",paste0("MS2_dataextraction_",Sys.Date(),"_Peaklists"),paste0("PC_",samplename,".csv")),row.names = F)
  
  n.ms2_spectraidentifier.11 <- n.ms2_spectraidentifier.10 %>% filter(PS_cont>=0.5) %>% select(Mass=precursorMZ,Time=retentionTime) %>% mutate(Mass=round(Mass,digits = 2),Time=round(x = Time,digits = 2),Name=c(paste0("PS_",seq_along((n.ms2_spectraidentifier.10 %>% filter(PS_cont>=0.5))$retentionTime))),Annotations="MS2derived")
  if (!dir.exists(paths = here("Output",paste0("MS2_dataextraction_",Sys.Date(),"_Peaklists")))){
    dir.create(path = here("Output",paste0("MS2_dataextraction_",Sys.Date(),"_Peaklists")))  }
  write.csv(n.ms2_spectraidentifier.11,file = here("Output",paste0("MS2_dataextraction_",Sys.Date(),"_Peaklists"),paste0("PS_",samplename,".csv")),row.names = F)
  getwd()     
}
