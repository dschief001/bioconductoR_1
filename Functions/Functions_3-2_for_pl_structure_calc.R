
# Functions used in pl_structure_calc.R -----------------------------------

###########################
#  FUNCTIONS DEFINITION   #
###########################

# a) Calculate masses of possible fragments FAs only ----------------------------------
fragments_theo <- function(carb, double) {
  Frag1_base <- 44.997654295    # getMass(getMolecule("HCO2")) # Define base m/z of fragmentation path 1
  CH2_base <- 14.01565007       # getMass(getMolecule("CH2")) # Define m/z difference for +1 carbon chain length
  H2_base <- 2.01565007         # getMass(getMolecule("H2")) # Define m/z difference for -1 double bond

  frag_theo <- c() # Initialize result object
  for(carcount in carb[1]:carb[2]){ # Open loop for side chain length from 20 to 50
    for(db in double[1]:double[2]){ # Open loop for 0 to 4 double bonds

      # Calculate the exact mass for specific fragment
      fragments <- data.frame(name = paste("FA ", carcount,":",db, sep=""),
                              carbons = carcount,
                              desat = db,
                              frag1 = Frag1_base + (carcount-1)*CH2_base  -db*H2_base)

      frag_theo <- rbind(frag_theo, fragments) # Write into result object
    }
  }
  fraglist <- cbind(melt(frag_theo, id.vars = c("name","carbons", "desat")), height = 0)
  fraglist <- fraglist[with(fraglist, order(fraglist$value)), ]
  fraglist <- as_tibble(fraglist)
  fraglist <- fraglist %>% unite(carbons,desat,sep = ":",col = "Name")
  return(fraglist)
}
# b) Calculate masses of possible fragments sn1 ether ----------------------------------
# 
# carb = range_carb
# double = range_double




fragments_theo_ether <- function(carb, double) {
  Frag1_base <- 31.018389735                    # getMass(getMolecule("H3CO")) # Define base m/z of fragmentation path 1
  CH2_base <-   14.01565007                     # getMass(getMolecule("CH2")) # Define m/z difference for +1 carbon chain length
  H2_base <-     2.01565007                     # getMass(getMolecule("H2")) # Define m/z difference for -1 double bond

  headgroup_PC_mitOH <-  225.076609712          # getMass(getMolecule("C7H16NO5P"))
  headgroup_PC_ohneOH <- 207.066045012          # getMass(getMolecule("C7H14NO4P"))

  headgroup_PE_mitOH <-  197.045309572          # getMass(getMolecule("C5H12NO5P"))
  headgroup_PE_ohneOH <- 179.034744872          # getMass(getMolecule("C5H10NO4P"))

  headgroupwHGL_PS_mitOH <- 154.003110395       # getMass(getMolecule("C3H7O5P"))
  headgroupwHGL_PS_ohneOH <- 135.992545695      # getMass(getMolecule("C3H5O4P"))

  frag_theo <- c() # Initialize result object
  frag_theo_midrange <- c() # Initialize result object

  for(carcount in carb[1]:carb[2]){ # Open loop for side chain length from 20 to 50
    for(db in double[1]:double[2]){ # Open loop for 0 to 4 double bonds

      # Calculate the exact mass for specific fragment
      sn1fragments <- data.frame(Name = paste("A ", carcount,":",db, sep=""),
                                 carbons = carcount,
                                 desat = db,
                                 mass = Frag1_base + (carcount-1)*CH2_base  -db*H2_base,
                                 Identifier = "sn1_Alcohol")

      mfragments1 <- data.frame(Name = paste("A ", carcount,":",db, sep=""),
                                carbons = carcount,
                                desat = db,
                                mass  = (Frag1_base + (carcount-1)*CH2_base  -db*H2_base)+headgroup_PC_mitOH,
                                Identifier = "sn1_alc+HG-PC_mitOH")
      mfragments2 <- data.frame(Name = paste("A ", carcount,":",db, sep=""),
                                carbons = carcount,
                                desat = db,
                                mass  = (Frag1_base + (carcount-1)*CH2_base  -db*H2_base)+headgroup_PC_ohneOH,
                                Identifier = "sn1_alc+HG-PC_ohneOH")
      mfragments3 <- data.frame(Name = paste("A ", carcount,":",db, sep=""),
                                carbons = carcount,
                                desat = db,
                                mass  = (Frag1_base + (carcount-1)*CH2_base  -db*H2_base)+headgroup_PE_mitOH,
                                Identifier = "sn1_alc+HG-PE_mitOH")
      mfragments4 <- data.frame(Name = paste("A ", carcount,":",db, sep=""),
                                carbons = carcount,
                                desat = db,
                                mass  = (Frag1_base + (carcount-1)*CH2_base  -db*H2_base)+headgroup_PE_ohneOH,
                                Identifier = "sn1_alc+HG-PE_ohneOH")
      mfragments5 <- data.frame(Name = paste("A ", carcount,":",db, sep=""),
                                carbons = carcount,
                                desat = db,
                                mass  = (Frag1_base + (carcount-1)*CH2_base  -db*H2_base)+headgroupwHGL_PS_mitOH,
                                Identifier = "sn1_alc-HG-PS_mitOH")
      mfragments6 <- data.frame(Name = paste("A ", carcount,":",db, sep=""),
                                carbons = carcount,
                                desat = db,
                                mass  = (Frag1_base + (carcount-1)*CH2_base  -db*H2_base)+headgroupwHGL_PS_ohneOH,
                                Identifier = "sn1_alc-HG-PS_ohneOH")
      #=> other midrange ether fragments add here
      frag_theo <- rbind(frag_theo, sn1fragments) # Write into result object
      mfragments <- rbind(mfragments1,mfragments2,mfragments3,mfragments4,mfragments5,mfragments6)
      frag_theo_midrange <- rbind(frag_theo_midrange, mfragments) #pass on to continuing function
    }
  }
  frag_theo <- rbind(frag_theo,frag_theo_midrange)
  colnames(frag_theo)
  fraglist <- cbind(melt(frag_theo, id.vars = c("Name","carbons", "desat","Identifier")), height = 0)
  fraglist <- fraglist[with(fraglist, order(fraglist$value)), ]
  fraglist <- as_tibble(fraglist)


  fraglist <- fraglist %>% mutate(Ceven= carbons %%2 ==0,Dbeven= desat %%2==0,possPlas=desat>=1)

  fraglist <- fraglist %>% unite(carbons,desat,sep = ":",col = "Name")
  return(fraglist)
}

# c) function to select the nearest peak for 3 calibration peaks ----------
select_nearest_peak <- function(df){
  rm(obs)
  masses <- df[,1]
  for (n in seq_along(cal_ref)) {
    if (any(min(abs(masses-cal_ref[n])) <= mass_tol_calibration)) {
      if (exists(x = "obs",where = -1)){
        obs <- add_row(.data = obs,cal_p=cal_ref[n],cal_o=masses[abs(masses-cal_ref[n]) <= mass_tol_calibration])
      } else {
        obs <- tibble(cal_p=cal_ref[n],cal_o=masses[abs(masses-cal_ref[n]) <= mass_tol_calibration])
      }
    }
  }
  if(!exists(x = "obs",where = -1)){
    obs <- tibble(cal_p=cal_ref[n],cal_o=masses[abs(masses-cal_ref[n]) <= mass_tol_calibration])
    obs <- NULL}
  return(obs)
}

# d) correct the masses with standardpeak ---------------------------------

correct_massdrift <- function(data) {
  colnames(data) <- c("V1","V2")
  data <- as_tibble(data)
  data <- add_column(data,correctedmass=predict(object = cal_fit,tibble(x=data$V1)))
  return(data)
}

# e) FA combinator with given FAs -----------------------------------------
# data <- n.ms2_spectraidentifier.2$data[[1]]
generate_possible_FA_layout <- function(data){
  
  #get carbon and DB info from Name, filter the generated FA distribution for the possible combinations only.
  
  fa <- c("10:0", "10:1",
          "12:0", "12:1",
          "13:0", "13:1",
          "14:0", "14:1", "14:2",
          "15:0", "15:1", "15:2",
          "16:0", "16:1", "16:2", "16:3", "16:4",
          "17:0", "17:1", "17:2", "17:3", "17:4",
          "18:0", "18:1", "18:2", "18:3", "18:4", "18:5", "18:6",
          "19:0", "19:1", "19:2", "19:3", "19:4", "19:5", "19:6",
          "20:0", "20:1", "20:2", "20:3", "20:4", "20:5", "20:6",
          "21:0", "21:1", "21:2", "21:3", "21:4", "21:5", "21:6",
          "22:0", "22:1", "22:2", "22:3", "22:4", "22:5", "22:6", "22:7","22:8",
          "23:0", "23:1", "23:2", "23:3", "23:4", "23:5", "23:6", "23:7","23:8",
          "24:0", "24:1", "24:2", "24:3", "24:4", "24:5", "24:6", "24:7","24:8",
          "25:0", "25:1", "25:2", "25:3", "25:4", "25:5", "25:6", "25:7","25:8",
          "26:0", "26:1", "26:2", "26:3", "26:4", "26:5", "26:6", "26:7","26:8",
          "27:0", "27:1", "27:2", "27:3", "27:4", "27:5", "27:6", "27:7","27:8",
          "28:0", "28:1", "28:2", "28:3", "28:4", "28:5", "28:6", "28:7","28:8"
  )
  
  # Alternatively a combinatorial fatty acid space can be constructed. This will result in a much larger number of degrees of freedom and much longer calculation times
  #fa <- paste(arrange(expand.grid(a=seq(10,30),b=seq(0,12)),a)[,1], arrange(expand.grid(a=seq(10,30),b=seq(0,12)),a)[,2], sep = ".")
  
  # Here all possible dual combination of fatty acids are constructed
  fa.comb <- as_tibble(dplyr::arrange(expand.grid(a=fa,b=fa),a))
  
  
  if (str_sub(string = unique(data$Name),start = 3,end = 3)=="P") {
          
    fa.comb <- fa.comb %>% separate(col = a,into = c("aC","aDB"),sep = ":") %>% separate(col = b,into = c("bC","bDB"),sep = ":") %>% filter(aDB>=1) %>% 
              mutate(cC=(as.numeric(aC)+as.numeric(bC)),cDB=(as.numeric(aDB)+as.numeric(bDB))) %>% filter(cC== unique(data$realC),cDB== unique(data$realDB))
  } else {
    fa.comb <- fa.comb %>% separate(col = a,into = c("aC","aDB"),sep = ":") %>% separate(col = b,into = c("bC","bDB"),sep = ":") %>% 
              mutate(cC=(as.numeric(aC)+as.numeric(bC)),cDB=(as.numeric(aDB)+as.numeric(bDB))) %>% filter(cC== unique(data$realC),cDB== unique(data$realDB))
  }
  fa.comb <- fa.comb %>% unite(aC,aDB,col="FA1",sep=":") %>% unite(bC,bDB,col="FA2",sep=":") %>% select(-c(cC,cDB))
}

# f) filter the FA_combinations table -------------------------------------
# data1 <- n.ms2_spectraidentifier.2$data[[1]]

actual_filterfunction <- function(data1){
  changedata <- data1$FA_Frag_mass[[1]]
  filtervalues <-data1$FA_possible[[1]]$FA1 
  changedata%>% filter(Name %in% filtervalues)
}
# g) filter the FA_combinations ETHER table -------------------------------------
actual_filterfunction_ether <- function(data1){
  changedata <- data1$Ether_restmass[[1]]
  filtervalues <-data1$FA_possible[[1]]$FA1 
  changedata%>% filter(Name %in% filtervalues)
}

# h) Mass_annotation function com bine FA to possible combinations --------
# View(n.ms2_spectraidentifier.8use)
# i=166
# n.ms2_spectraidentifier.8use[165,]
# test <- for(i in seq_along(n.ms2_spectraidentifier.8use$data)) {
#   data <- n.ms2_spectraidentifier.8use$data[[i]]
#   data$MS2_data_masscorr
#   ms2_annotation(data)
#   print(i)
#   flush.console()
# }



ms2_annotation <- function(data){
  
  FA_fragment_masses_saveinto <- data$FA_Frag_massN[1][[1]]
  #select midrange fragments...
  tester_class <- str_split_fixed(string = unique(data$Name),pattern = ":",n = 2)[,1]
  tester_class <- str_sub(tester_class,end = 2)
  M_fragment_masses_saveinto <- data$M_Frag_massN[1][[1]]
  selM_fragment_masses_saveinto <- M_fragment_masses_saveinto %>% filter(str_detect(string = Identifier,pattern = tester_class))
  
  measured_masses_gethere <- data$MS2_data_masscorr[1][[1]]
  # working with the measured MS2 peaks here:
  # therefore calculate the tester sums by massrange here:
  tot_FA_signal <- measured_masses_gethere %>% filter((min(data$FA_Frag_mass[[1]]$value)-1) < correctedmass ) %>% filter(correctedmass < (max(data$FA_Frag_mass[[1]]$value)+1))
  filtercrit <- data$Ether_restmass[[1]] %>% filter(!str_detect(string = Identifier,pattern = "sn1_Alcohol"))
  tot_M_signal <- measured_masses_gethere %>% filter((min(filtercrit$value)-1) < correctedmass ) %>% filter(correctedmass < (max(filtercrit$value)+1))
  
  
  mass_tol
  mass_tol2
  
  for (mu in seq_along(FA_fragment_masses_saveinto$Name)) {
    if(nrow(measured_masses_gethere)>0) {
      if(abs(measured_masses_gethere[which(abs(FA_fragment_masses_saveinto$value[mu] - measured_masses_gethere$correctedmass) == min(abs(FA_fragment_masses_saveinto$value[mu] - measured_masses_gethere$correctedmass))),3]-FA_fragment_masses_saveinto$value[mu]) < mass_tol){
        # Is the fragmented precursor inside the mass tolerance?
        FA_fragment_masses_saveinto[mu,5] <- measured_masses_gethere[which(abs(FA_fragment_masses_saveinto$value[mu] - measured_masses_gethere$correctedmass) == min(abs(FA_fragment_masses_saveinto$value[mu] - measured_masses_gethere$correctedmass))),2] # If yes, exract the number of the respective fragment spectrum
      }
    }
  } 
  for (mu in seq_along(selM_fragment_masses_saveinto$Name)) {
    if(nrow(measured_masses_gethere)>0) {
      if(abs(measured_masses_gethere[which(abs(selM_fragment_masses_saveinto$value[mu] - measured_masses_gethere$correctedmass) == min(abs(selM_fragment_masses_saveinto$value[mu] - measured_masses_gethere$correctedmass))),3]-selM_fragment_masses_saveinto$value[mu]) < mass_tol2){
        # Is the fragmented precursor inside the mass tolerance?
        selM_fragment_masses_saveinto[mu,5] <- measured_masses_gethere[which(abs(selM_fragment_masses_saveinto$value[mu] - measured_masses_gethere$correctedmass) == min(abs(selM_fragment_masses_saveinto$value[mu] - measured_masses_gethere$correctedmass))),2] # If yes, exract the number of the respective fragment spectrum
      }
    }
  } 
  FA_Indicator <- sum(FA_fragment_masses_saveinto$height)/sum(tot_FA_signal$V2)
  M_Indicator <- sum(selM_fragment_masses_saveinto$height)/sum(tot_M_signal$V2)
  # negativ beispiel
  # Peparatatory guess Lipidclass from data
  identificator <- measured_masses_gethere %>% filter(!(round(data$precursorMZ,digits= 0) == round(measured_masses_gethere$V1,digits = 0)) & V1 > lower_threshold_Identification)
  (data$precursorMZ - 60)
  
  PC_cont=(identificator %>% filter(mass_tol_identification >= abs(identificator$V1-(data$precursorMZ - 60))))$V2/sum(identificator$V2)
  PS_cont=(identificator %>% filter(mass_tol_identification >= abs(identificator$V1-(data$precursorMZ - 87))))$V2/sum(identificator$V2)
  PC_count=(identificator %>% filter(mass_tol_identification >= abs(identificator$V1-(data$precursorMZ - 60))))$V2
  PS_count=(identificator %>% filter(mass_tol_identification >= abs(identificator$V1-(data$precursorMZ - 87))))$V2
  
  Find_PC_PS <- tibble(Name=data$Name,Sample=samplename,precursorMZ=data$precursorMZ,
                       PC_cont = ifelse(length(PC_cont)>0,PC_cont,NA),
                       PS_cont = ifelse(length(PS_cont)>0,PS_cont,NA),
                       PC_count =ifelse(length(PC_count)>0,PC_count,NA),
                       PS_count =ifelse(length(PS_count)>0,PS_count,NA))
  # sum(identificator$V2)
  # 
  # if(abs(data$precursorMZ - (identificator %>% filter(V2 ==max(V2)))$V1) > 87-mass_tol_identification &
  #   abs(data$precursorMZ - (identificator %>% filter(V2 ==max(V2)))$V1) < 87+mass_tol_identification){
  #     class.specific <- (identificator %>% filter(V2 == max(V2)) %>% summarise(Value=sum(V2)))$Value
  #     a <- tibble(Class="PS",Value = c(class.specific/(sum(identificator$V2))))
  # }else{
  #   a <- 0
  #   }
  # if(abs(data$precursorMZ - (identificator %>% filter(V2 ==max(V2)))$V1) > 60-mass_tol_identification &
  #   abs(data$precursorMZ - (identificator %>% filter(V2 ==max(V2)))$V1) < 60+mass_tol_identification){
  #     class.specific <- (identificator %>% filter(V2 == max(V2)) %>% summarise(Value=sum(V2)))$Value
  #     a <- tibble(Class="PS",Value = c(class.specific/(sum(identificator$V2))))
  # }else{
  #   a <- 0
  #   }
  # 
  # 
  # if(87 == round(data$precursorMZ - (identificator %>% filter(V2 ==max(V2)))$V1,digits = 0)) {
  #   print("is PS")
  #   class.specific <- identificator %>% filter(V2 == max(V2)) %>% summarise(sum(V2))
  #   a <- tibble(Class="PS",Value = c(class.specific$`sum(V2)`/(sum(identificator$V2))))
  #   identificator2 <- identificator %>% filter(V2 != max(V2))
  #   b <- ()
  #   if(60 == round(data$precursorMZ - (identificator2 %>% filter(V2 == max(V2)))$V1,digits = 0)){
  #     class.specific <- identificator2 %>% filter(V2 == max(V2)) %>% summarise(sum(V2))
  #     b <- tibble(Class="PC",Value = c(class.specific$`sum(V2)`/(sum(identificator$V2))))
  #   }
  #   bind_rows(a,b)
  #   
  #   identificator3 <- identificator2 %>% filter(V2 != max(V2))
  #   round(data$precursorMZ - (identificator3 %>% filter(V2 == max(V2)))$V1,digits = 0)
  #   
  #   
  #   
  #   
  #   }
  # 
  # if(60 == round(data$precursorMZ - (identificator %>% filter(V2 ==max(V2)))$V1,digits = 0)) {
  #   print("is PC")
  #   identificator <- identificator %>% filter(V2 != max(V2))
  #   round(data$precursorMZ - (identificator %>% filter(V2 == max(V2)))$V1,digits = 0)
  #   
  #   }
  
  return(list(FA_fragment_masses_saveinto,selM_fragment_masses_saveinto,FA_Indicator,M_Indicator,Find_PC_PS))
}

# i) Diagnostic evaluation of the data: -----------------------------------

diagnostic_data <- function(data) {
  df <- tibble(FA_diag=data[[3]],M_diag = data[[4]])
}

remove_diagdata <- function(data){
  discard(.x = data,.p = c(F,F,T,T,T))
}

# j) Clean up data --------------------------------------------------------
cleanup_FA_possible <- function(x){
  str_sort  <- function(x) {x %>% sort() %>% str_c(collapse = "-")}
  TEMP.2 <- apply(X = as.data.frame(x), MARGIN = 1, FUN=function(x) str_sort(x)) %>% 
    unique() %>% str_split(pattern = "-") %>% unlist() %>% matrix(ncol = 2, byrow = T)
}

# k) here a feedback predictor function: ---------------------------------
# data <- n.ms2_spectraidentifier.9
Feedback_plots <- function(data){
  Feedback.1<- data %>% ungroup() %>%  transmute(Name1_,seqNum,Feedback=map(Found_FAs,.f = diagnostic_data)) %>% unnest(cols = c(Feedback))
  Feedback.2 <- Feedback.1 %>% mutate(FA_diag2=if_else(is.na(FA_diag)==T,true = -1,false = FA_diag),M_diag2=if_else(is.na(M_diag)==T,true = -1,false = M_diag))
  Feedback.3.1 <- arrange(.data = Feedback.2,seqNum) %>% mutate(seqNum=as.factor(seqNum))
  d1 <- ggplot(data = Feedback.3.1,mapping = aes(x = seqNum,y = FA_diag2))+
    geom_col()+
    geom_hline(mapping = aes(yintercept = 1,col="red"))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90,vjust=0.5), legend.position = "none")+
    labs(title = "FA fragments coverage",caption = "normal [0-1], -1 = no data covered")
  d4 <- ggplot(data = Feedback.3.1,mapping = aes(x = Name1_,y = seqNum,fill=Name1_))+
    geom_tile()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90,vjust=0.5),legend.position = "none")+
    labs(title = "FA MS2 coverage total")
  
  Indicator.1 <- summary(Feedback.3.1$seqNum)[summary(Feedback.3.1$seqNum)>=2]
  Feedback.3.2 <- Feedback.3.1 %>% filter(seqNum %in% names(Indicator.1))
  d5 <- ggplot(data = Feedback.3.2,mapping = aes(x = Name1_,y = seqNum,fill=Name1_))+
    geom_tile()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90,vjust=0.5,size = 12),legend.position = "none")+
    labs(title = "FA MS2 Overlap")
  
  d2 <- ggplot(data = Feedback.3.1,mapping = aes(x = seqNum,y = M_diag2))+
    geom_col()+
    geom_hline(mapping = aes(yintercept = 1,col="red"))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90,vjust=0.5), legend.position = "none")+
    labs(title = "midrange fragments coverage",caption = "normal [0-1], -1 = no data covered")
  
  
  if (!dir.exists(paths = here("output",paste("MS2_dataextraction_",Date_today)))) {
    dir.create(path = here("output",paste("MS2_dataextraction_",Date_today)))}
  if (!dir.exists(paths = here("output",paste("MS2_dataextraction_",Date_today),"2.Overlap"))) {
    dir.create(path = here("output",paste("MS2_dataextraction_",Date_today),"2.Overlap"))}
  if (!dir.exists(paths = here("output",paste("MS2_dataextraction_",Date_today),"1.Diagnostic"))) {
    dir.create(path = here("output",paste("MS2_dataextraction_",Date_today),"1.Diagnostic"))}
  
  d3 <- arrangeGrob(grobs = list(d1,d2),nrow = 2,ncol = 1)
  d6 <- arrangeGrob(grobs = list(d5,d4),nrow = 2,ncol = 1)
  ggsave(plot = d3,filename = paste0(samplename,"_MS2_dataextraction_Diagnostic.pdf"),path = here("Output",paste("MS2_dataextraction_",Date_today),"1.Diagnostic"),width = 30,height = 10)
}

# l) Shuffle the single FA data into the combinatoric prepared emp --------
## testing setup for following function ##
# rm("Name","a","b","c","e")
# LAST_FA_profiles_intens <- list()
# uo=1
# for (uo in seq_len(length.out = nrow(n.ms2_spectraidentifier.12))) {
#   print(uo)
#   Name=n.ms2_spectraidentifier.12$Name[[uo]]
#   a=n.ms2_spectraidentifier.12$FA_possible[[uo]]
#   b=n.ms2_spectraidentifier.12$FA_possible2[[uo]]
#   c= n.ms2_spectraidentifier.12$Found_FAs[[uo]]
#   e=n.ms2_spectraidentifier.12$totIonCurrent[[uo]]
#   LAST_FA_profiles_intens[[uo]] <- sort_FAS_into_emtytable_multiplyPrecIntens(Name,a,b,c,e)
#   rm("Name","a","b","c","e")
# }
# rm("Name","a","b","c","e")
# Name=n.ms2_spectraidentifier.12$Name[[4]]
# a=n.ms2_spectraidentifier.12$FA_possible[[4]]
# b=n.ms2_spectraidentifier.12$FA_possible2[[4]]
# c= n.ms2_spectraidentifier.12$Found_FAs[[4]]
# e=n.ms2_spectraidentifier.12$totIonCurrent[[4]]
# sort_FAS_into_emtytable_multiplyPrecIntens(Name,a,b,c,e)

sort_FAS_into_emtytable_multiplyPrecIntens <- function(Name,a,b,c,e) {
  
  Subclass <- str_split_fixed(string = Name,pattern = ":",n = 2)[,1]
  composition <- str_split_fixed(string = Name,pattern = ":",n = 2)[,2]
  
  # filter decarboxylation products
  FA_meas <- c[[1]]
  FA_meas.1 <- FA_meas %>% arrange(desc(height))
  FA_meas.2 <- FA_meas.1
  z=1
  for (z in 1:5) {
    if (any(abs(round(FA_meas.2$value,digits = 0) - round(FA_meas.2$value[z]-44,digits = 0))<=0.3)) {
      FA_meas.2 <- FA_meas.2 %>% slice(-which(abs(round(FA_meas.2$value,digits = 0) - round(FA_meas.2$value[z]-44,digits = 0))<=0.3))
    }
  }
  FA_meas <- FA_meas.2
  
  #distinguish between acyl and alkyl/alk-1-enyl linkage by MS1 assignment
  if (str_detect(string = Subclass,pattern = "(?<=[A-Z]{2})[OpeP]{1}")) {
    # ether lipid processing MS2
    Assigned <- str_extract(string = composition,pattern = "[^:]*:[^:]*")
    poss_FAs <- unique(c(a[,1], a[,2]))
    b <- tibble(sn2_FA=poss_FAs,height=0)
    #FA assignment to sn2 position only
    for(q in seq_along(FA_meas$Name)){
      if(any(str_detect(string = poss_FAs,pattern = FA_meas$Name[q]))){
        b$height[b$sn2_FA == as.character(FA_meas$Name[q])] <- as.numeric(FA_meas$height[q])
        #anteilig aufgeteilt durch div by occurence
      }else{
        b$height[b$sn2_FA == as.character(FA_meas$Name[q])] <- 0
      }
    }
    ret_data <- list()
    ret_data[[1]] <- b
    names(ret_data) <- "FA_sn2_data"
    #calculated the sn2_FA_residue contributions
    b1 <- arrange(.data = b,desc(sn2_FA))
    b1$height <- 0
    colnames(b1) <- c("sn1_FA","height2")
    b2 <- b1
    
    # #confirm with markus
    # test <- tibble(sn2=b$sn2_FA,b2$sn2_FA)
    
    M_meas <- c[[2]]
    M_meas.2 <- M_meas %>% filter(Name %in% poss_FAs)
    M_meas.3 <- M_meas.2
    
    # add ether treatment:
    # need to filter possible found peaks only!
    #maybe a bit hard filtering but will test like this
    # x <- (M_meas.2 %>% group_by(Name) %>% nest())$data[[19]]
    M_meas.2 <- M_meas.2 %>% group_by(Name) %>% nest() %>% mutate(data2=map(.x = data,.f = function(x){
      x <- x %>% arrange(Identifier)
      x$height[1]>x$height[2]*1.5 #& (all(x$height>0))
    })) %>% unnest(cols = c("data","data2")) %>% filter(data2==T) %>% select(-data2)
    
    if (nrow(M_meas.2)==0) {
      M_meas.2 <- M_meas.3 %>% mutate(height = 0)
    }
    
    for (q in seq_along(unique(M_meas.2$Identifier))) {
      assign(paste0(if_else(condition = str_detect(string = unique(M_meas.2$Identifier)[q],pattern = "ohne"),true = "ohne",false = "mit"),"_OH"),M_meas.2 %>% filter(Identifier==unique(M_meas.2$Identifier)[q]))
    }
    #without OH
    if ("ohne_OH" %in% ls()) {
      for(q in seq_along(ohne_OH$Name)){
        
        if(any(str_detect(string = b1$sn1_FA ,pattern = ohne_OH$Name[q]))){
          b1$height2[b1$sn1_FA == as.character(ohne_OH$Name[q])] <- as.numeric(ohne_OH$height[q])
          #anteilig aufgeteilt durch div by occurence
        }else{
          b1$height2[b1$sn1_FA == as.character(ohne_OH$Name[q])] <- 0
        }
      }
    }else{stop(paste(Subclass,composition,"error in ohne_OH treatment,
                     Class seems to need midrangefragments, seem undefined"))}
    # withOH
    for(q in seq_along(mit_OH$Name)){
      
      if(any(str_detect(string = b2$sn1_FA ,pattern = mit_OH$Name[q]))){
        b2$height2[b2$sn1_FA == as.character(mit_OH$Name[q])] <- as.numeric(mit_OH$height[q])
        #anteilig aufgeteilt durch div by occurence
      }else{
        b2$height2[b2$sn1_FA == as.character(mit_OH$Name[q])] <- 0
      }
    }
    
    #order and sort together possible combinations of fatty acids at sn1 and sn2 position.
    # hier eine liste und mit map arbeiten??
    ohne_OH <- cbind(arrange(b,sn2_FA),arrange(b1, desc(sn1_FA))) %>% select(sn1_FA,sn2_FA,sn1_height=height2,sn2_height=height)
    mit_OH <- cbind(arrange(b,sn2_FA),arrange(b2, desc(sn1_FA))) %>% select(sn1_FA,sn2_FA,sn1_height=height2,sn2_height=height)
    
    ret_data[[2]] <- ohne_OH
    ret_data[[3]] <- mit_OH
    names(ret_data)[2:3] <- c("ohne_OH","mit_OH")
    
    m_ohne_OH <- ohne_OH %>% gather(sn1_FA,sn2_FA,key = "position",value = "FA") %>% transmute(position,FA,height=(sn1_height+sn2_height))
    m_mit_OH <- mit_OH %>% gather(sn1_FA,sn2_FA,key = "position",value = "FA") %>% transmute(position,FA,height=(sn1_height+sn2_height))
    
    
    ret_data[[4]] <- m_ohne_OH
    ret_data[[5]] <- m_mit_OH
    names(ret_data)[4:5] <- c("m_ohne_OH","m_mit_OH")
    # 1. FA an sn2 position weighted only
    
    tester <- sum(ret_data[[1]]$height)!=0
    if(tester){
      temp <- ret_data[[1]] %>% mutate(Profile0_1= height/sum(height))
      ret_data[[6]] <- temp %>% mutate(Intensity_prof=Profile0_1*e)
    } else{
      temp <- ret_data[[1]]  %>% mutate(Profile0_1= height)
      ret_data[[6]] <- temp %>% mutate(Intensity_prof=Profile0_1*e) 
    }
    # 2. weighting of averaged compositions with/without OH
    tester <- sum(ret_data[[4]]$height)!=0
    if(tester){
      temp <- ret_data[[4]] %>% group_by(position) %>% mutate(Profile0_1= height/sum(height))
      ret_data[[7]] <- temp %>% mutate(Intensity_prof=Profile0_1*e)
    } else{
      temp <- ret_data[[4]] %>% group_by(position) %>% mutate(Profile0_1= height)
      ret_data[[7]] <- temp %>% mutate(Intensity_prof=Profile0_1*e) 
    }
    tester <- sum(ret_data[[5]]$height)!=0 
    if(tester){
      temp <- ret_data[[5]] %>% group_by(position) %>% mutate(Profile0_1= height/sum(height))
      ret_data[[8]] <- temp %>% mutate(Intensity_prof=Profile0_1*e)
    } else{
      temp <- ret_data[[5]] %>% group_by(position) %>% mutate(Profile0_1= height)
      ret_data[[8]] <- temp %>% mutate(Intensity_prof=Profile0_1*e) 
    }
    
    tester <- sum(ret_data[[4]]$height)!=0 & sum(ret_data[[5]]$height) != 0
    if(tester){
      
      temp.1 <- ret_data[[4]]
      temp.1$sample <- "ohneOH"
      temp.2 <- ret_data[[5]]
      temp.2$sample <- "mitOH"
      temp.3 <- bind_rows(temp.1,temp.2)
      temp.4 <- temp.3 %>% group_by(position,FA) %>% summarise(height=sum(height))
      temp <- temp.4 %>% group_by(position) %>% mutate(Profile0_1= height/sum(height))
      ret_data[[9]] <- temp %>% mutate(Intensity_prof=Profile0_1*e)
    } else{
      temp <- ret_data[[5]] %>% group_by(position) %>% mutate(Profile0_1= height)
      ret_data[[9]] <- temp %>% mutate(Intensity_prof=Profile0_1*e) 
    }
    
    names(ret_data)[6:9] <- c("w_FA_sn2_data","w_ohne_OH","w_mit_OH","mean_W_FA")
    return(ret_data)
    # wieder rowsums und mitteln.... um auszugleichen die unterschiedliche intensität... trotzdem zuordnung behalten.
    #   - normalisiert auch das Verhältnis
  }
  else{
    # normal PL processing als done already
    FA_meas <- c[[1]]
    for(q in seq_along(FA_meas$Name)){
      
      if(any(str_detect(string = a,pattern = FA_meas$Name[q]))){
        b[a == as.character(FA_meas$Name[q])] <- as.numeric(FA_meas$height[q])/table(a == as.character(FA_meas$Name[q]))[2]
        #anteilig aufgeteilt durch div by occurence
      }else{
        b[a == as.character(FA_meas$Name[q])] <- 0
      }
    }
    d <- tibble(FA=c(a[,1], a[,2]),Value=c(rowSums(b)/2, rowSums(b)/2))
    d <- d %>% group_by(FA) %>% summarise(Abundance = sum(Value)) 
    tester <- sum(d$Abundance)!=0
    if(tester){
      d <- d %>% mutate(Profile0_1= Abundance/sum(Abundance))
      d %>% mutate(Intensity_prof=Profile0_1*e) 
    } else{
      d <- d%>% mutate(Profile0_1= Abundance)
      d %>% mutate(Intensity_prof=Profile0_1*e) 
    }
  }
}

# m) plotting function MS2 data -------------------------------------------

# data=n.ms2_spectraidentifier.13$LAST_FA_profiles_intens[[1]]
plot_FA_results <- function(data,seqNum,Name){
  Subclass <- str_split_fixed(string = Name,pattern = ":",n = 2)[,1]
  Assigned <- str_extract(string = str_split_fixed(string = Name,pattern = ":",n = 2)[,2],pattern = "[^:]*:[^:]*")
  if (!is.data.frame(data)) {
    # 1. ether treatment
    
    ###################### FA raw data
    e.1 <- ggplot(data = data$FA_sn2_data,mapping = aes(x = sn2_FA,y = height))+
      geom_col(width = .8)+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5))+
      labs(title = paste("seqNum",seqNum,"--",Subclass,Assigned,"*sn2* FA composition ether lipid"),caption = "this resembles the overall FA distribution at sn1 and sn2 position, not taking the exact binding site into account",
           y="Counts [au]",x="Fatty acid species")
    ################## ether lipid raw data
    e.2.1 <- ggplot(data = data$ohne_OH,mapping = aes(x = sn2_FA,y = sn2_height))+
      geom_col(width = .8,fill="dodgerblue3")+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "none")+
      labs(title = paste("seqNum",seqNum,"--",Subclass,Assigned,"raw [e/p] no OH"),
           y="Counts [au]",x="sn2 Acyl residues")
    
    e.2.2 <- ggplot(data = data$ohne_OH,mapping = aes(x = sn1_FA,y = sn1_height))+
      geom_col(width = .8,fill="darkorange2")+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "none")+
      labs(title = paste("seqNum",seqNum,"--",Subclass,Assigned,"raw [e/p] no OH"),
           y="Counts [au]",x="sn1 Alkyl residues")
    
    e.2.3 <- ggplot(data = data$mit_OH,mapping = aes(x = sn2_FA,y = sn2_height))+
      geom_col(width = .8,fill="dodgerblue3")+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "none")+
      labs(title = paste("seqNum",seqNum,"--",Subclass,Assigned,"raw [e/p] with OH"),caption = "residues raw data resolved according to binding site",
           y="Counts [au]",x="sn2 Acyl residues")
    
    e.2.4 <- ggplot(data = data$mit_OH,mapping = aes(x = sn1_FA,y = sn1_height))+
      geom_col(width = .8,fill="darkorange2")+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "none")+
      labs(title = paste("seqNum",seqNum,"--",Subclass,Assigned,"raw [e/p] with OH"),
           y="Counts [au]",x="sn1 Alkyl residues")
    e.2 <-arrangeGrob(e.2.1,e.2.2,e.2.3,e.2.4,top = "raw data ether lipids")
    ###################### ether lipid averaged over binding site data
    e.3.1 <- ggplot(data = data$m_ohne_OH,mapping = aes(x = FA,y = height,fill=position))+
      geom_col(width = .8)+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "right")+
      labs(title = " no OH",caption = "residues averaged signals per combination colour corresponds with binding site",
           y="Counts [au]",x="residues")+
      scale_fill_manual(values = c("#E69F00", "#56B4E9"),labels=c("sn1 alkyl/alkenyl","sn2 acyl"))
    e.3.2 <- ggplot(data = data$m_mit_OH,mapping = aes(x = FA,y = height,fill=position))+
      geom_col(width = .8)+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "right")+
      labs(title = "with OH",caption = "residues averaged signals per combination colour corresponds with binding site",
           y="Counts [au]",x="residues")+
      scale_fill_manual(values = c("#B22222", "#228B22"),labels=c("sn1 alkyl/alkenyl","sn2 acyl"))
    
    e.3 <- arrangeGrob(e.3.1,e.3.2, top=paste("seqNum",seqNum,"--",Subclass,Assigned,"**averaged** [e/p]","averaged ether lipid profiles"))
    ############### abundance plots of averaged data
    e.4.1 <- ggplot(data = data$w_ohne_OH,mapping = aes(x = FA,y = Profile0_1 ,fill=position))+
      geom_col(width = .8)+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "right")+
      labs(title = " no OH",
           y="Abundance",x="residues")+
      scale_fill_manual(values = c("#DA70D6", "#00CED1"),labels=c("sn1 alkyl/alkenyl","sn2 acyl"))
    
    e.4.2 <- ggplot(data = data$w_mit_OH,mapping = aes(x = FA,y = Profile0_1 ,fill=position))+
      geom_col(width = .8)+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "right")+
      labs(title = "with OH",caption = "residues averaged signals per combination colour corresponds with binding site",
           y="Abundance",x="residues")+
      scale_fill_manual(values = c("#FF1493", "#32CD32"),labels=c("sn1 alkyl/alkenyl","sn2 acyl"))
    
    e.4 <- arrangeGrob(e.4.1,e.4.2, top=paste("seqNum",seqNum,"--",Subclass,Assigned,"**relative** [e/p]","averaged ether lipid profiles"))
    ################################## total plot averaged on count data  
    e.5 <- ggplot(data = data$mean_W_FA,mapping = aes(x = FA,y = Profile0_1 ,fill=position))+
      geom_col(width = .8)+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "right")+
      labs(title = paste("seqNum",seqNum,"--",Subclass,Assigned),caption = "total residues averaged signals per combination colour corresponds with binding site",
           y="Abundance",x="residues")+
      scale_fill_manual(values = c("#FF1493", "#32CD32"),labels=c("sn1 alkyl/alkenyl","sn2 acyl"))
    
    Ms2_plots <- list(e.1,e.2,e.3,e.4,e.5)
    return(Ms2_plots)
    
  }else{
    # ester phospholipids treatment
    e.1 <- ggplot(data = data,mapping = aes(x = FA,y = Profile0_1))+
      geom_col(width = .8)+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5))+
      labs(title = paste("seqNum",seqNum,"--",Subclass,Assigned,"FA composition averaged over phospholipid"),caption = "this resembles the overall FA distribution at sn1 and sn2 position, not taking the exact binding site into account",
           y="Abundance [%]",x="Fatty acid species")
    Ms2_plots <- list(e.1)
    return(Ms2_plots)
  }
  
}


# n) final plotting function writes plots to pdf separated by plot --------

plotting_function <- function(data) {
  # 1. create folder for plots
  if (!dir.exists(paths = here("output",paste("MS2_dataextraction_",Date_today)))){
    dir.create(path = here("output",paste("MS2_dataextraction_",Date_today)))  }
  if (!dir.exists(paths = here("output",paste("MS2_dataextraction_",Date_today),samplename))) {
    dir.create(path = here("output",paste("MS2_dataextraction_",Date_today),samplename))}
  # 2. separate Name for plotting
  Subclass <- str_split_fixed(string = data$Name,pattern = ":",n = 2)[,1]
  composition <- str_split_fixed(string = data$Name,pattern = ":",n = 2)[,2]
  Assigned <- str_extract(string = composition,pattern = "[^:]*:[^:]*")
  #distinguish between acyl and alkyl/alk-1-enyl linkage by MS1 assignment
  
  # 3. ether lipid processing MS2
  {
    Ether_plots <- data %>% filter(str_detect(string = Name,pattern = "(?<=[A-Z]{2})[OPpe]{1}"))
    temp.1 <- list()
    temp.2 <- list()
    temp.3 <- list()
    temp.4 <- list()
    temp.5 <- list()
    for (z in seq_along(Ether_plots$Name)) {
      temp.1[[z]] <- Ether_plots$plotlist[[z]][[1]]
      temp.2[[z]] <- Ether_plots$plotlist[[z]][[2]]
      temp.3[[z]] <- Ether_plots$plotlist[[z]][[3]]
      temp.4[[z]] <- Ether_plots$plotlist[[z]][[4]]
      temp.5[[z]] <- Ether_plots$plotlist[[z]][[5]]
    }  
    p1 <- marrangeGrob(grobs = temp.1,ncol = 1,nrow = 2,bottom = quote(paste("page", g, "of",npages)))
    ggsave(plot = p1,filename = paste0("ETHER-1_MS2_sn2_FA_raw_data.pdf"),path = here("Output",paste("MS2_dataextraction_",Date_today),samplename),width = 8.3,height = 11.7)                                                           
    
    p1 <- marrangeGrob(grobs = temp.2,ncol = 1,nrow = 2,bottom = quote(paste("page", g, "of",npages)))
    ggsave(plot = p1,filename = paste0("ETHER-2_MS2_both_raw_data.pdf"),path = here("Output",paste("MS2_dataextraction_",Date_today),samplename),width = 8.3,height = 11.7)                                                           
    
    p1 <- marrangeGrob(grobs = temp.3,ncol = 1,nrow = 2,bottom = quote(paste("page", g, "of",npages)))
    ggsave(plot = p1,filename = paste0("ETHER-3_MS2_both_averaged.pdf"),path = here("Output",paste("MS2_dataextraction_",Date_today),samplename),width = 8.3,height = 11.7)                                                           
    
    p1 <- marrangeGrob(grobs = temp.4,ncol = 1,nrow = 2,bottom = quote(paste("page", g, "of",npages)))
    ggsave(plot = p1,filename = paste0("ETHER-4_MS2_percentualcontributions_averaged.pdf"),path = here("Output",paste("MS2_dataextraction_",Date_today),samplename),width = 8.3,height = 11.7)                                                           
    
    p1 <- marrangeGrob(grobs = temp.5,ncol = 1,nrow = 4,bottom = quote(paste("page", g, "of",npages)))
    ggsave(plot = p1,filename = paste0("ETHER-5_MS2_TOTAL_percentual-averaged.pdf"),path = here("Output",paste("MS2_dataextraction_",Date_today),samplename),width = 8.3,height = 11.7)                                                           
    
  }
  # 4. ester lipid processing for plotting
  Ester_plots <- data %>% filter(str_detect(string = Name,pattern = "(?<=^[A-Z]{2})[OpeP]{1}",negate = T))
  temp.1 <- list()
  temp.1 <- foreach(z = seq_along(Ester_plots$Name), .export = c("Ester_plots")) %dopar% {
    Ester_plots$plotlist[[z]][[1]]
  }  
  p1 <- marrangeGrob(grobs = temp.1,ncol = 1,nrow = 2,bottom = quote(paste("page", g, "of",npages)))
  ggsave(plot = p1,filename = paste0("ESTER-1_MS2_sn2_FA_raw_data.pdf"),path = here("Output",paste("MS2_dataextraction_",Date_today),samplename),width = 8.3,height = 11.7)                                                           
}


# o) preprare final intensity data ----------------------------------------
filter_average_MS2_data <- function(intens){
  if (is.data.frame(intens)) {
    intens
  }else{
    intens$mean_W_FA
  }
}

# p) weighting of the data ------------------------------------------------
# Name = n.ms2_F.2$Name[[1]]
# data = n.ms2_F.2$data[[1]]

averaging_over_spectra <-  function(Name,data){
  rm(Subclass)
  Subclass <- str_split_fixed(string = Name,pattern = ":",n = 2)[,1]
  
  if (str_detect(string = Subclass,pattern = "(?<=^[A-Z]{2})[OpeP]{1}")) {
    temp <- bind_rows(data$LAST_FA_profiles_intens2)
    temp.2 <- temp %>% group_by(position, FA)%>% summarise(FA_sum=sum(Intensity_prof))
    temp.2 %>% ungroup()%>% mutate(FA_norm=(FA_sum/sum(FA_sum)))
    
  }else{
    temp <- bind_rows(data$LAST_FA_profiles_intens2)
    temp.2 <- temp %>% group_by(FA)%>% summarise(FA_sum=sum(Intensity_prof))
    temp.2 %>% ungroup() %>% mutate(FA_norm=(FA_sum/sum(FA_sum)))
  }
}


# q) plotting function 2, averaged per sample and peak --------------------
# data <- n.ms2_F.4

plotting_function2 <- function(data) {
  # 1. create folder for plots
  if (!dir.exists(paths = here("output",paste("MS2_dataextraction_",Date_today)))){
    dir.create(path = here("output",paste("MS2_dataextraction_",Date_today)))  }
  if (!dir.exists(paths = here("output",paste("MS2_dataextraction_",Date_today),samplename))) {
    dir.create(path = here("output",paste("MS2_dataextraction_",Date_today),samplename))}
  # 2. separate Name for plotting
  data.1 <- data %>% group_by(Name) %>% nest()
  Subclass <- str_split_fixed(string = data$Name,pattern = ":",n = 2)[,1]
  
  Assigned <- str_extract(string = str_split_fixed(string = data$Name,pattern = ":",n = 2)[,2],pattern = "[^:]*:[^:]*")
  #distinguish between acyl and alkyl/alk-1-enyl linkage by MS1 assignment
  
  # 3. ether lipid processing MS2
  {
    Ether_plots <- data.1 %>% filter(str_detect(string = Name,pattern = "(?<=^[A-Z]{2})[OpeP]{1}"))
    
    temp.1 <- list()
    library("doFuture")
    registerDoFuture()
    # plan(multiprocess)
    plan(sequential)
    
    temp.1 <- foreach(z = seq_along(Ether_plots$Name), .export = c("Ether_plots")) %dopar% {
      Subclass <- str_split_fixed(string = Ether_plots$Name[[z]],pattern = ":",n = 2)[,1]
      Assigned <- str_extract(string = str_split_fixed(string = Ether_plots$Name[[z]],pattern = ":",n = 2)[,2],pattern = "[^:]*:[^:]*")
      
      x <- ggplot(data = Ether_plots$data[[z]],mapping = aes(x = FA,y = FA_norm ,fill=position))+
        geom_col(width = .8)+
        theme_bw()+
        theme(axis.text.x = element_text(angle=90,vjust=0.5),legend.position = "right")+
        labs(title = paste("FA composition averaged over phospholipid",Subclass,Assigned),caption = "total residues averaged signals per combination colour corresponds with binding site",
             y="Abundance",x="residues")+
        scale_fill_manual(values = c("#FF1493", "#32CD32"),labels=c("sn1 alkyl/alkenyl","sn2 acyl"))
      return(x)
      }  
    p1 <- marrangeGrob(grobs = temp.1,ncol = 1,nrow = 2,bottom = quote(paste("page", g, "of",npages)))
    ggsave(plot = p1,filename = paste0("ETHER-lipid_composition_",samplename,".pdf"),path = here("Output",paste("MS2_dataextraction_",Date_today),samplename),width = 8.3,height = 11.7)                                                           
  }
  # 4. ester lipid processing for plotting
  temp.1 <- list()
  Ester_plots <- data.1 %>% filter(str_detect(string = Name,pattern = "(?<=^[A-Z]{2})[OpeP]{1}",negate = T))
  temp.1 <- foreach(z = seq_along(Ester_plots$Name), .export = c("Ester_plots")) %dopar% {
    Subclass <- str_split_fixed(string = Ester_plots$Name[[z]],pattern = ":",n = 2)[,1]
    Assigned <- str_extract(string = str_split_fixed(string = Ester_plots$Name[[z]],pattern = ":",n = 2)[,2],pattern = "[^:]*:[^:]*")
    
    x <- ggplot(data = Ester_plots$data[[z]],mapping = aes(x = FA,y = FA_norm))+
      geom_col(width = .8)+
      theme_bw()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5))+
      labs(title = paste("FA composition averaged over phospholipid",Subclass,Assigned),caption = "this resembles the overall FA distribution at sn1 and sn2 position, not taking the exact binding site into account",
           y="Abundance [%]",x="Fatty acid species")
    return(x)
    }  
  p1 <- marrangeGrob(grobs = temp.1,ncol = 1,nrow = 2,bottom = quote(paste("page", g, "of",npages)))
  ggsave(plot = p1,filename = paste0("ESTER-composition_",samplename,".pdf"),path = here("Output",paste("MS2_dataextraction_",Date_today),samplename),width = 8.3,height = 11.7)                                                           
}


# o) get contribution data out --------------------------------------------
# x=n.ms2_spectraidentifier.9

get_contribution_data <- function(x){
  x <- x %>% mutate(Contrib=map(Found_FAs,.f = function(x) {x[[5]]}))
  bind_rows(x$Contrib)
}

# p) Contribution data saved into overlap folder of the correct ti --------
Contribution_plot <- function(Contribution){
  if (!dir.exists(paths = here("output",paste("MS2_dataextraction_",Date_today),"2.Overlap"))) {
    dir.create(path = here("output",paste("MS2_dataextraction_",Date_today),"2.Overlap"))}
  temp <- Contribution %>% gather(PC_cont,PS_cont,key = "Class",value = "Value")
  d1 <- ggplot(data = temp,mapping = aes(x = Name,y = Value,fill=Class))+
    geom_col()+
    theme_bw()+
    facet_wrap(~Sample,ncol = 1)+
    theme(axis.text.x = element_text(angle=90,size=10,vjust=0.5))  
  ggsave(plot = d1,filename = paste0(samplename,"_MS2_dataextraction_Contribution.pdf"),path = here("Output",paste("MS2_dataextraction_",Date_today),"2.Overlap"),width = 8,height = 5)
}
