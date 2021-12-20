# Wed Jan 23 10:06:21 2019 ------------------------------


# Functions for new integrated data processing ----------------------------

# Feedback funtion for Batchidentifier ------------------------------------
Feedback_Batch_identifier <- function(input) {
  if(str_detect(string = Batch_identifier,pattern = "_")) {
    stop(msg="Batch_identifier contains not allowed Zeichen! [_]
       1. check Batch_identifier input
       2. restart processing")
  }else{
    print("Well done")}
}


# Data readout function ---------------------------------------------------
# data_input <- mzmine.integration$File[1]

Reslist_Rtlist_readout <- function(data_input){
  Filename <- str_split(string = data_input,pattern = "\\\\(?=[^\\\\]*$)",simplify = T)[,2]
  data <- tryCatch(expr = read.csv(file = data_input,header = T),error=function(e){print(paste(Filename," <= data not readable"))})
  
  
  
  temp <- tryCatch(tibble(Name=as.character(data[,grep(pattern = ".ID.$",x = colnames(data))]),
                          Area=as.numeric(data[,grep(pattern = "Peak.area$",x = colnames(data))]),
                          tStart=data[,grep(pattern = "RT.start$",x = colnames(data))],
                          tEnd=data[,grep(pattern = "RT.end$",x = colnames(data))],
                          Mass=data[,grep(pattern = "row.m.z{1}",x = colnames(data))],
                          Comment=as.character(data[,grep(pattern = "row.comment",x = colnames(data))])),
                   
                   error=function(e) {paste("Error in the readout of data in sample",Filename)})
  tryCatch(expr = {
    temp %>% mutate(Batch=Batch_identifier,Class=str_split(string = Filename,pattern = "_",simplify = T,n=2)[1],
                    Sample= str_remove(string = str_split(string = Filename,pattern = "_",simplify = T,n=2)[2],
                                       pattern = "\\.csv"))
  },error=function(e) {paste("Error in the readout of data in sample and adding Identifier columns in",Filename)})
}


# quantification function -------------------------------------------------
# data <- mzmine.quantification.2$data[[2]]
# Class <- mzmine.quantification.2$Class[[2]]

quantification <- function(Class,data){
  if (all(str_detect(string = data$Name,pattern = "^[PS]{1}[ACEGSIM]{1}[Ope]*(?=\\:.)"))& length(unique(str_sub(string = data$Name,end = 2)))==1) {
    PL.class <- unique(unique(str_sub(string = data$Name,end = 2)))
  } else {
    ## b) via drawer item names 
    print(paste("Check",Class,"for wrongly assigned peak or more than one Class / wrong class?!!",
                unique(str_extract(string = data$Name, pattern = "^[PS]{1}[ACEGSIM]{1}[Ope]*(?=\\:.)"))))
  }
  if (PL.class == "SM") {
    PL.class <- "PC"
  }
  if (PL.class == "PI") {
    PL.class <- "PE"
  }
  if (all(PL.class != "SM"&!is.na(PL.class)&PL.class != "PI")) {
    #get the corresponding cal.table values:
    quad_i.koeff <- cal.table %>% filter(Class == PL.class)
    
    #Mit quadratischer Lösungsformel umrechnen in Stoffmengen [mol]
    data <- data %>% filter(Area > 0) %>% mutate(logArea=log10(Area),Checker=is.infinite(logArea))
    data <- data %>% mutate(quant_mol=(10^(
      ((-quad_i.koeff$b+sqrt(quad_i.koeff$b^2-4*quad_i.koeff$a*(quad_i.koeff$c-logArea)))/(2*quad_i.koeff$a)))*Dilution_factor))
    #HIER BERÜCKSICHTIGT, dass nur 10µL der gesamt extrahierten Menge initiiert wurden. dh zehnfache Menge extrahiert pro gemessenem Proteinwert, daher *10
  }
}


# Extract cumulative ISTD signals per class and tissue --------------------

ISTD_readout_classwise <- function(data){
    PL.class <- unique(unique(str_sub(string = data$Name,end = 2)))
    criteria <- all(!(PL.class %in% c("SM","NA" ,"PI")))
    if (criteria) {
      ISTD.data <- data %>% filter(str_detect(string = data$Name,pattern = paste0(PL.class,":28:0"))) %>% group_by(Sample) %>% summarise(ISTD=sum(quant_mol))
    }
}

# Normalisation function to ISTD extraction efficiency --------------------
# x = mzmine.quantification.3$quantified[[1]]
# y = mzmine.quantification.3$ISTD[[1]]

normalisation_ISTD <- function(x,y){
  left_join(x = x,y = y,by="Sample") %>% mutate(normed_mol=(quant_mol/ISTD)*ISTD_amount)
}

# Normalisation function to protein content -------------------------------

normalisation_protein_content <- function(x){
  x %>% mutate(normed_mol_mg_protein=(normed_mol/Protein_values_mg))
}

########## Setup general parameters -------------

Dilution_factor <- 10
ISTD_amount <- 5.00e-10 # mol in resolubilised lipid extract 
pl.pattern <- "P[ACEGIS]{1}[ep]*|SM"
Quantification_pattern <- "^[PS]{1}[ACEGSIM]{1}[Ope]*(?=\\:.)"


