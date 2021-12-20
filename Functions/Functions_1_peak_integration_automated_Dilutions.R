# additional Functions ---------------------------------------------------------------

##' Peakfinder
# does a peak shape optimisation in time domain.

peakfinder <- function(chrom, cutoff){
  #cutoff=0.1
  midpoint <- round((length(chrom)+1)/2, digits = 0)
  # plot(chrom)
  delta <- ifelse(test = midpoint< 15,yes = midpoint-2,no = 15)


  chrom_b <- chrom-max(max(chrom[(midpoint-delta):(midpoint+delta)],na.rm = T))*cutoff >= 0
  chrom_b <- c(F, chrom_b, F)

  index<- midpoint+1

  while (chrom_b[index] == T) {
    index <- index - 1
  }
  low <- index

  index<- round((length(chrom)+1)/2, digits = 0)
  while (chrom_b[index] == T) {
    index <- index + 1
  }
  high <- index-2


  if(high-low < 40){ # Change 40 to adjust for min peak length (40=0.37 min)
    diff <- 40-(high-low)
    low <- low - floor(diff/2)
    if(low<1){low<-1}
    high <- high + ceiling(diff/2)
    if(high>length(chrom)){high<-length(chrom)}

  }

  return(c(low, high))
}


#' Peakfinder 2
#'
#' This function allows to integrate peaks starting from the peak maximum and then extening the "peakwindow" timewise to the front and end of the feature until the intensitiy value drops below a specified threshold called `cutoff`.
#' The input parameter `chrom` should contain the extracted sprectra already filtered by massrange and possible retentiontime.
#'
#' @param chrom
#' @param cutoff
#'
#' @return
#' @export
#'
#' @examples
peakfinder2 <- function(chrom, cutoff){
  #cutoff=0.1
  chromx <- chrom
  # as.numeric((chromx %>% filter(Peak_max==T))$Spectrum)
  chrom <-  chrom$Intensity
  midpoint <- (chromx %>% rowid_to_column(var = "ID") %>% filter(Peak_max==T))$ID

  chrom_b <- chrom-(chromx %>% filter(Peak_max==T))$Intensity*cutoff >= 0
  length(chrom_b)
  # test
  chrom_b <- c(F, chrom_b, F)

  index <- midpoint
  while (chrom_b[index] == T) {
    index <- index - 1
  }
  low <- index

  index <- midpoint
  while (chrom_b[index] == T) {
    index <- index + 1
  }
  high <- index-2

  min_peak_length <- 10 # just guessing!

  if(high-low < min_peak_length){ # Change 40 to adjust for min peak length (40 = 0.37 min)
    diff <- min_peak_length-(high-low)
    low <- low - floor(diff/2)
    if(low<1){low<-1}
    high <- high + ceiling(diff/2)
    if(high>length(chrom)){high<-length(chrom)}

  }

  return(c(low, high))
}

#
#' Peakfinder based on slopes
#'
#' Allows to detect feature start and end by change of the slope. sensitivity can be changed by variation of the `multiplier` value.
#'
#' @param chrom
#' @param multiplier
#'
#' @return
#' @export
#'
#' @examples
peakfinder3 <- function(chrom, multiplier){
  #cutoff=0.1

  chromx <- chrom
  # as.numeric((chromx %>% filter(Peak_max==T))$Spectrum)
  chrom <-  chrom$Intensity
  midpoint <- (chromx %>% rowid_to_column(var = "ID") %>% filter(Peak_max==T))$ID

  chromx[(midpoint-10):midpoint,4]
  vector <- c()
  nrow(chromx)-midpoint
  for (n in 1:(midpoint-10)) {
    temp <- chromx[(midpoint-10-n):(midpoint-n),]
    x =temp$rt_times
    y = temp$Intensity
    model <- lm(y~x)
    vector[n]= model$coefficients[2]
  }
  vector2 = vector >= max(vector)*multiplier
  vector2 <- c( vector2, F)
  vector
  index <- 1
  while (vector2[index] == T) {
    index <- index + 1
  }
  low <- midpoint - index-1

  for (n in 1:(nrow(chromx)-midpoint-1)) {
    temp <- chromx[(midpoint+n):(midpoint+10+n),]
    x =temp$rt_times
    y = temp$Intensity
    model <- lm(y~x)
    vector[n]= model$coefficients[2]

  }

  vector2 = vector <= max(vector,na.rm = T)*multiplier
  vector2 <- c( vector2, F)
  vector
  index <- 1
  while (vector2[index] == T) {
    index <- index + 1
  }
  high <- midpoint + index-1
  return(c(low,high))
}
## maybe here a minimum peaklength needs to be defined also?
# min_peak_length <- 10 # just guessing!
#
# if(high-low < min_peak_length){ # Change 40 to adjust for min peak length (40 = 0.37 min)
#   diff <- min_peak_length-(high-low)
#   low <- low - floor(diff/2)
#   if(low<1){low<-1}
#   high <- high + ceiling(diff/2)
#   if(high>length(chrom)){high<-length(chrom)}
#
# }
#
# return(c(low, high))
# }


# Peakfinder 4
# include shoulders:
#
peakfinder4 <- function(chrom, multiplier,multiplier2){
  #cutoff=0.1
  multiplier = 0.01
  multiplier2 = 0.06
  chromx <- chrom3
  chromx <- chrom
  # as.numeric((chromx %>% filter(Peak_max==T))$Spectrum)
  chrom <-  chrom3$Intensity
  chrom <-  chrom$Intensity

  midpoint <- (chromx %>% rowid_to_column(var = "ID") %>% filter(Peak_max==T))$ID

  # chromx[(midpoint-10):midpoint,4]
  vector <- c()
  nrow(chromx)-midpoint
  for (n in 1:(midpoint-10)) {
    temp <- chromx[(midpoint-10-n):(midpoint-n),]
    x =temp$rt_times
    y = temp$Intensity
    model <- lm(y~x)
    vector[n]= model$coefficients[2]
  }
  vector2 = vector >= max(vector)*multiplier
  vector2 <- c( vector2, F)
  # vector
  index <- 1
  while (vector2[index] == T) {
    index <- index + 1
  }
  low <- midpoint - index-1
  vector <- c()
  for (n in 1:(nrow(chromx)-midpoint-1)) {
    temp <- chromx[(midpoint+n):(midpoint+10+n),]
    x =temp$rt_times
    y = temp$Intensity
    model <- lm(y~x)
    vector[n]= model$coefficients[2]

  }

  vector2 = vector <= max(vector,na.rm = T)*multiplier
  vector2 <- c( vector2, F)
  vector
  index <- 1
  while (vector2[index] == T) {
    index <- index + 1
  }
  high <- midpoint + index-1


  if ((chromx$Intensity[low]-chromx$Intensity[high]) >= 800) {
    print("go to lower")
    n=0
    y = chromx$Intensity[low]
    while (abs(y-chromx$Intensity[high])>= 1000 & n+10>=1) {
      n = n-1
      temp <- chromx[(low+n):(low-10+n),]
      y = temp$Intensity
      y = mean(y)
    }
    low <- low+n
  }
  # try to again find step decent:
  if ((chromx$Intensity[high]-chromx$Intensity[low]) >= 800) {
    ## linear k estimation loop, not implemented to find peak end here...
    # for (n in 1:(nrow(chromx)-high-1)) {
    #   temp <- chromx[(high+n):(high+10+n),]
    #   x =temp$rt_times
    #   y = temp$Intensity
    #   model <- lm(y~x)
    #   vector[n]= model$coefficients[2]
    # }
    #
    n=0
    y = chromx$Intensity[high]
    while (abs(y-chromx$Intensity[low])>= 1000 & n+10<=(nrow(chromx)-high-1)) {
      n = n+1
      temp <- chromx[(high+n):(high+10+n),]
      y = temp$Intensity
      y = mean(y)
    }
    high <- n+high

    chromx$rt_times[high]

    vector <- c()

    for (n in seq(from=high-20,to = high-10)) {
      temp <- chromx[(n):(n+10),]
      x =temp$rt_times
      y = temp$Intensity
      model <- lm(y~x)
      vector[n]= model$coefficients[2]
    }
    vector <- vector[!is.na(vector)]
    # barplot(vector)
    vector2 = abs(vector) >= max(abs(vector),na.rm = T)*multiplier2
    vector2 <- c( vector2, F)
    index <- 1
    while (vector2[index] == T) {
      index <- index + 1
    }

    high <- high-20+index-1
  }
  #
  #
  # }
  #
  #   chromx$Intensity[low]
  # chromx$Intensity[high]
  #
  #
  return(c(low,high))
}

# based on gaussian peak shape from midpoint of the peak.
peakfinder5 <- function(chrom) {
  #cutoff=0.1
  # chromx <- chrom3
  chromx <- chrom
  # as.numeric((chromx %>% filter(Peak_max==T))$Spectrum)
  # chrom <-  chrom3$Intensity
  chrom <-  chrom$Intensity

  FWHM_height <- (chromx %>% rowid_to_column(var = "ID") %>% filter(Peak_max==T))$Intensity*0.5

  temp <- (chromx %>% mutate(dev=sign(FWHM_height-Intensity)))$dev
  midpoint <- (chromx %>% rowid_to_column(var = "ID") %>% filter(Peak_max==T))$ID

  cpv <- temp[midpoint]
  index <- 0
  while (temp[midpoint-index] == cpv && midpoint-index > 1) {
    index <- index + 1
  }
  low <- midpoint - index

  index <- 0
  while (temp[midpoint+index] == cpv && index+midpoint < length(temp)) {
    index <- index + 1
  }
  high <- midpoint + index

  sigma=(chromx$rt_times[high]-chromx$rt_times[low])/2.38
  y <- qnorm(p = c(0.05,0.95), mean=chromx$rt_times[midpoint], sd=sigma)

  y  # would already contain the correct gaussian times! for 90% of te peak integrated.
  # added this part to ensure a minimum peaklength
  if(abs(y[2]-y[1])<=0.12){
    y[1] <- y[1]-0.05
    y[2] <- y[2]+0.05
  }
  return(c(which.min(abs(chromx$rt_times-y[1])),which.min(abs(chromx$rt_times-y[2]))))
}
## Peakfinder 6: enlarged gaussian shaped time window... by lowering FWHM estimation value...

peakfinder6 <- function(chrom) {
  #cutoff=0.1
  # chromx <- chrom3
  chromx <- chrom
  # as.numeric((chromx %>% filter(Peak_max==T))$Spectrum)
  # chrom <-  chrom3$Intensity
  chrom <-  chrom$Intensity

  FWHM_height <- (chromx %>% rowid_to_column(var = "ID") %>% filter(Peak_max==T))$Intensity*0.48

  temp <- (chromx %>% mutate(dev=sign(FWHM_height-Intensity)))$dev
  midpoint <- (chromx %>% rowid_to_column(var = "ID") %>% filter(Peak_max==T))$ID

  cpv <- temp[midpoint]
  index <- 0
  while (temp[midpoint-index] == cpv && midpoint-index > 1) {
    index <- index + 1
  }
  low <- midpoint - index

  index <- 0
  while (temp[midpoint+index] == cpv && index+midpoint < length(temp)) {
    index <- index + 1
  }
  high <- midpoint + index

  sigma=(chromx$rt_times[high]-chromx$rt_times[low])/2.38
  y <- qnorm(p = c(0.05,0.95), mean=chromx$rt_times[midpoint], sd=sigma)

  y  # would already contain the correct gaussian times! for 90% of te peak integrated.
  # added this part to ensure a minimum peaklength
  if(abs(y[2]-y[1])<=0.12){
    y[1] <- y[1]-0.05
    y[2] <- y[2]+0.05
  }
  return(c(which.min(abs(chromx$rt_times-y[1])),which.min(abs(chromx$rt_times-y[2]))))
}
## Peakfinder 7: enlarged gaussian shaped time window... by lowering FWHM estimation value...

peakfinder7 <- function(chrom) {
  #cutoff=0.1
  # chromx <- chrom3
  chromx <- chrom
  # as.numeric((chromx %>% filter(Peak_max==T))$Spectrum)
  # chrom <-  chrom3$Intensity
  chrom <-  chrom$Intensity

  FWHM_height <- (chromx %>% rowid_to_column(var = "ID") %>% filter(Peak_max==T))$Intensity*0.5

  temp <- (chromx %>% mutate(dev=sign(FWHM_height-Intensity)))$dev
  midpoint <- (chromx %>% rowid_to_column(var = "ID") %>% filter(Peak_max==T))$ID

  cpv <- temp[midpoint]
  index <- 0
  while (temp[midpoint-index] == cpv && midpoint-index > 1) {
    index <- index + 1
  }
  low <- midpoint - index

  index <- 0
  while (temp[midpoint+index] == cpv && index+midpoint < length(temp)) {
    index <- index + 1
  }
  high <- midpoint + index

  sigma=(chromx$rt_times[high]-chromx$rt_times[low])/2.38
  y <- qnorm(p = c(0.01,0.99), mean=chromx$rt_times[midpoint], sd=sigma)

  y  # would already contain the correct gaussian times! for 90% of te peak integrated.
  # added this part to ensure a minimum peaklength
  if(abs(y[2]-y[1])<=0.12){
    y[1] <- y[1]-0.05
    y[2] <- y[2]+0.05
  }
  return(c(which.min(abs(chromx$rt_times-y[1])),which.min(abs(chromx$rt_times-y[2]))))
}



# checking the slope variant:
# library(foreach)
# listi <- foreach (i =seq(from=5,to = 50,by = 5),.combine = "rbind") %do% {
#   for (n in 1:30) {
#   temp <- chromx[(midpoint+n):(midpoint+i+n),]
#     x =temp$rt_times
#     y = temp$Intensity
#   model <- lm(y~x)
#   vector[n]= model$coefficients[2]
#   }
# # tibble(1:30,vector)
# return(tibble(x=1:30,slope=vector,section=i))
# vector <- numeric(length = 30)
# }
#
# ggplot(listi,mapping = aes(x = x,y = slope,group=section,color=section))+geom_line()

# use the here got indixes to calculate the regressions!




# Peak selection:

# selects the highest peak:
peakmax <- function(x, mass){
  x[which.min(abs(x[,1] - mass)), 2]
}

# Selects the highest Peak at a local m/z maximum closest to the given m/z value
# if no peak is found the corresponding time is removed from the list.

# x <- object$Spectra[[1]] # %>% filter(m_z > 857 & m_z < 859) %>% view()
# x <- spectra[[100]]

peakmax2 <- function(x, mass){
  # plot(x)
  # plot(x,xlim = c(1429.95,1430.2))
  y <- x %>% filter(abs(m_z-mass)<0.4)
  #
  # plot(y)
  # abline(v = mass)
  y <- y %>% filter(if_else(lag(Intensity) < Intensity & lead(Intensity) < Intensity, TRUE, FALSE)) %>% mutate(Xc=abs(m_z-mass)) %>% filter(Xc==min(Xc)) %>% select(-Xc)
  if (nrow(y)==0) {x %>% filter(abs(m_z-mass)==min(abs(m_z-mass)))} else{y}
}

# Get relevant region from raw data
mzfinder2 <- function(x, peak_mz_range){
  x <- x[x[,1] >= peak_mz_range[1] & x[,1] <= peak_mz_range[2], ]
  tibble(m_z=x[ ,1],Intensity=x[ ,2])
}


### calibrate function:
# takes an input mass and time and searches in the proximities for the highest Intensity (mz direction the local maximum nearest to the given mass).
# Returns the TIC table per peak.


# x <- rtstdlist[3,]
# mz=mz
# ms1_header=ms1_header


# Functions 2 for Peak integration ----------------------------------------

# Extracts the calibration peak info, with TIC in time domain.
calibrate <- function(x, mz, ms1_header){
  # Delta_t <- 2.0 #min
  Delta_t <- 0.8 #min
  lipid = x[1]
  rt = as.numeric(x[2])
  mass = as.numeric(x[3])

  rt_min <-  ms1_header[which.min(abs(ms1_header$retentionTime - (rt*60 - Delta_t*60))), 1]
  rt_max <-  ms1_header[which.min(abs(ms1_header$retentionTime - (rt*60 + Delta_t*60))), 1]
  rt_spectra <- ms1_header$seqNum[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]
  rt_times <- ms1_header$retentionTime[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]

  spectra <- peaks(mz, rt_spectra)

  # chrom <- unlist(lapply(X = spectra, FUN = function(x) peakmax(x, mass=mass)))

  object <- tibble(Spectra=spectra) %>% mutate(Spectra=map(Spectra,.f = function(x) {tibble(m_z=x[,1],Intensity=x[,2])}))
  object <- object %>% transmute(chrom_list=map(Spectra,.f = function(x) peakmax2(x, mass=mass)))
  chrom2 <- object %>% unnest(cols = c(chrom_list))

  # a modification of the time list, if no local maximum was found, then the entry is removed!
  if(any(map(.x = object$chrom_list,.f = nrow)==0)){
    rt_times2 <- rt_times[-c(which(unlist(map(.x = object$chrom_list,.f = nrow)==0)))]
  }else{
    rt_times2 <- rt_times
  }

  # plot(chrom)
  # plot(chrom2$Intensity)
  # plot(chrom2$m_z)
  # hist(chrom2$m_z,main = paste(mass))
  # return(rt_times[which.max(chrom)]/60)
  chrom3 <- chrom2 %>% mutate(Peak_max=Intensity==max(Intensity),rt_times=rt_times2/60,lipid=lipid)
  chrom3[which(round(chrom3$rt_times,digits = 5)==12.49372),]
  # return(chrom3$rt_times[which.max(chrom3$Intensity)])
  return(chrom3)
}




#
# i=
# for (i in 1:173) {
#   x <- lipidlist[i,]
#   peakextrakt2(rt=as.numeric(x[4]),
#                mass=as.numeric(x[3]),
#                mz=mz,
#                ms1_header=ms1_header,
#                ms2_header= ms2_header)
#   print(i)
# }
#
# x <- lipidlist[1,]
# x <- lipidlist[1,]
# rt=as.numeric(x[4])
# mass=as.numeric(x[3])
# mz=mz
# ms1_header=ms1_header

#' Peakfinder function
#' 
#' used to integrate the inputed peaks to provide an area readout per indicated feature in the feature list.
#'
#' @param rt 
#' @param mass 
#' @param mz 
#' @param ms1_header 
#' @param ms2_header 
#'
#' @return
#' @export
#'
#' @examples
peakextrakt2 <- function(rt, mass, mz, ms1_header,ms2_header){

  Delta_t <- 0.15 #min 60*0.3
  # print(mass)
  # rt=as.numeric(lipidlist[1, 4])
  # mass=as.numeric(lipidlist[1, 3])

  rt_min <-  ms1_header[which.min(abs(ms1_header$retentionTime - (rt*60 - Delta_t*60))), 1]
  rt_max <-  ms1_header[which.min(abs(ms1_header$retentionTime - (rt*60 + Delta_t*60))), 1]
  rt_spectra <- ms1_header$seqNum[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]
  rt_times <- ms1_header$retentionTime[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]

  spectra <- peaks(mz, rt_spectra)
  object <- tibble(Spectra=spectra) %>% mutate(Spectra=map(Spectra,.f = function(x) {tibble(m_z=x[,1],Intensity=x[,2])}))
  object <- object %>% transmute(chrom_list=map(Spectra,.f = function(x) peakmax2(x, mass=mass)))
  chrom2 <- object %>% unnest(cols = c(chrom_list))

  # a modification of the time list, if no local maximum was found, then the entry is removed!
  if(any(map(.x = object$chrom_list,.f = nrow)==0)){
    rt_times2 <- rt_times[-c(which(unlist(map(.x = object$chrom_list,.f = nrow)==0)))]
  }else{
    rt_times2 <- rt_times
  }

  chrom3 <- chrom2 %>% mutate(Peak_max=Intensity==max(Intensity),rt_times=rt_times2/60)
  # plot(chrom3$Intensity)
  # ggplot(data = chrom3,mapping = aes(x = rt_times,y = m_z,fill=Intensity))+
  #     geom_tile()
  ggplot(data = chrom3,mapping = aes(x = rt_times,y = Intensity))+
    geom_line()+geom_point()
  # Returns a matrix where each row represents one peak found. The first column gives the height, the second the position/index where the maximum is reached, the third and forth the indices of where the peak begins and ends --- in the sense of where the pattern starts and ends.

  rt2 <- (chrom3 %>% filter(Peak_max==T))$rt_times # rt time in min
  ################# adapt to shift window to maximum peak found.
  # adaption of the peakwindow to found maximum, allowed once, if the rt shifts are below a tolerance?
  decider <- rt-rt2 # Retentiontime (min)
  Delta_t2 <- 0.3
  rt_min <-  ms1_header[which.min(abs(ms1_header$retentionTime - (rt2*60 - Delta_t2*60))), 1]
  rt_max <-  ms1_header[which.min(abs(ms1_header$retentionTime - (rt2*60 + Delta_t2*60))), 1]
  rt_spectra <- ms1_header$seqNum[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]
  rt_times <- ms1_header$retentionTime[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]

  spectra <- peaks(mz, rt_spectra)
  object <- tibble(Spectra=spectra) %>% mutate(Spectra=map(Spectra,.f = function(x) {tibble(m_z=x[,1],Intensity=x[,2])}))
  object <- object %>% transmute(chrom_list=map(Spectra,.f = function(x) peakmax2(x, mass=mass)))
  chrom2 <- object %>% unnest(cols = c(chrom_list))

  # a modification of the time list, if no local maximum was found, then the entry is removed!
  if(any(map(.x = object$chrom_list,.f = nrow)==0)){
    rt_times2 <- rt_times[-c(which(unlist(map(.x = object$chrom_list,.f = nrow)==0)))]
  }else{
    rt_times2 <- rt_times
  }
  ## either take the max peak as starting point, or the found peak in the first rt time adjustment step (second line)
  # chrom3 <- chrom2 %>% mutate(Peak_max=Intensity==max(Intensity),rt_times=rt_times2/60)
  chrom3 <- chrom2 %>% mutate(rt_times=rt_times2/60,Peak_max= rt2==rt_times)
  #
  # chrom3 %>% filter(Peak_max==T)
  # ONLY CHANGE here! ---------------
  # plot(chrom3$Intensity)
  # rt_box <- peakfinder2(chrom = chrom3,cutoff = 0.1) # newer smaller rt window may result
  rt_box <- peakfinder3(chrom = chrom3,multiplier = 0.01) # newer smaller rt window may result [use 0.01 as standard for multiplier]
  # rt_box <- c(47,64)
  # rt_box <- peakfinder4(chrom = chrom3,multiplier = 0.01,multiplier2 = .01) # window always minimum of 40
  # rt_box <- peakfinder5(chrom = chrom3) # window always minimum of 40
  # rt_box <- peakfinder6(chrom = chrom3) # window always minimum of 40
  # rt_box <- peakfinder7(chrom = chrom3) # window always minimum of 40
  rt_box
  chrom3c <- chrom3[rt_box[1]:rt_box[2],]
  # plot(chrom3c$Intensity)

  # look at these plots for how good the peakintegration works...

  # ggplot(data = chrom3c,mapping = aes(x = rt_times,y = m_z,fill=Intensity))+
  #   geom_tile()
  # ggplot(data = chrom3,mapping = aes(x = rt_times,y = Intensity))+
  #   geom_line()+geom_point()+
  #   geom_area(data = chrom3c,mapping = aes(x = rt_times,y = Intensity,fill="red",alpha=0.6))+theme(legend.position = "none")
  # # print("optimisation done")
  peak_spectra_range <- rt_spectra[rt_box[1]:rt_box[2]]
  peak_rt_range <- rt_times2[rt_box[1]:rt_box[2]]
  # Extract optimal mz range

  # spectrum <- spectra[[which.max(chrom)]]
  spectrum <- spectra[[which.max(chrom3$Intensity)]]
  #plot(spectrum, type = "l")
  section <- (spectrum[spectrum[,1] > (chrom3 %>% filter(Peak_max==T))$m_z-0.3 & spectrum[,1] < (chrom3 %>% filter(Peak_max==T))$m_z+0.3, ])
  #plot(section, type = "l")
  # Add a more sophisticated method here

  rt_spectra2 <- ms2_header[ms2_header$seqNum >= min(peak_spectra_range) & ms2_header$seqNum <= max(peak_spectra_range),]
  rt_spectra2 <- rt_spectra2 %>% filter(abs(chrom3 %>% filter(Peak_max==T) %>% pull(m_z)-rt_spectra2$precursorMZ)<0.3) %>% select(seqNum,precursorMZ)

  MS2 <- peaks(mz, rt_spectra2 %>% pull(seqNum))
  names(MS2) <- rt_spectra2 %>% pull(precursorMZ)
  # print("MS2 done")

  peak_mz_range <- c(min(section[,1]), max(section[,1]))

  subset1 <- map(.x = peaks(mz, peak_spectra_range),.f = function(x) mzfinder2(x, peak_mz_range=peak_mz_range))
  subset2 <- map(.x = peaks(mz, peak_spectra_range),.f = function(x) mzfinder2(x, peak_mz_range=peak_mz_range-1.0078))
  subset3 <- map(.x = peaks(mz, peak_spectra_range),.f = function(x) mzfinder2(x, peak_mz_range=peak_mz_range+1.0078))
  peakmatrix1 <- bind_rows(subset1,.id = "Spectrum")
  peakmatrix2 <- bind_rows(subset2,.id = "Spectrum")
  peakmatrix3 <- bind_rows(subset3,.id = "Spectrum")

  return(chrom3 %>% filter(Peak_max==T) %>% select(D_m_z=m_z,D_rt_times = rt_times,D_rt_times=rt_times) %>%
           mutate(C_mz_low=peak_mz_range[1]
                  ,C_mz_high=peak_mz_range[2]
                  ,C_rt_Start = rt_times2[rt_box[1]]/60
                  ,C_rt_End=rt_times2[rt_box[2]]/60
                  ,B_Area=sum(peakmatrix1$Intensity)
                  ,B_M_1=sum(peakmatrix3$Intensity)
                  ,B_Blank=sum(peakmatrix2$Intensity)) %>%
           select(B_Area,B_Blank,everything()) %>% mutate(
             Peakshape = list(list(chrom3,chrom3c)),C_small_rt_dev=decider
             ,F_Fragmentation=list(MS2)
             )
  )
}

# used for second script!! check with JK
peak_reinintegrate <- function(rt_Start, rt_End,mass,  mz, ms1_header){
  Delta_t <- 0.15 #min 60*0.3
  # print(mass)
  # rt=as.numeric(lipidlist[1, 4])
  # mass=as.numeric(lipidlist[1, 3])

  rt_min <-  ms1_header[which.min(abs(ms1_header$retentionTime - (rt_Start*60))), 1]
  rt_max <-  ms1_header[which.min(abs(ms1_header$retentionTime - (rt_End*60))), 1]
  rt_spectra <- ms1_header$seqNum[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]
  rt_times <- ms1_header$retentionTime[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]

  spectra <- peaks(mz, rt_spectra)
  object <- tibble(Spectra=spectra) %>% mutate(Spectra=map(Spectra,.f = function(x) {tibble(m_z=x[,1],Intensity=x[,2])}))
  object <- object %>% transmute(chrom_list=map(Spectra,.f = function(x) peakmax2(x, mass=mass)))
  chrom2 <- object %>% unnest(cols = c(chrom_list))

  # a modification of the time list, if no local maximum was found, then the entry is removed!
  if(any(map(.x = object$chrom_list,.f = nrow)==0)){
    rt_times2 <- rt_times[-c(which(unlist(map(.x = object$chrom_list,.f = nrow)==0)))]
  }else{
    rt_times2 <- rt_times
  }

  chrom3 <- chrom2 %>% mutate(Peak_max=Intensity==max(Intensity),rt_times=rt_times2/60)

  # ggplot(data = chrom3,mapping = aes(x = rt_times,y = Intensity))+
  #   geom_line()+geom_point()
  # Returns a matrix where each row represents one peak found. The first column gives the height, the second the position/index where the maximum is reached, the third and forth the indices of where the peak begins and ends --- in the sense of where the pattern starts and ends.
  # integrate whole area as input defines...

  rt_box <- c(1,nrow(chrom3))
  chrom3c <- chrom3[rt_box[1]:rt_box[2],]
  # plot(chrom3c$Intensity)

  # look at these plots for how good the peakintegration works...

  # ggplot(data = chrom3c,mapping = aes(x = rt_times,y = m_z,fill=Intensity))+
  #   geom_tile()
  # ggplot(data = chrom3,mapping = aes(x = rt_times,y = Intensity))+
  #   geom_line()+geom_point()+
  #   geom_area(data = chrom3c,mapping = aes(x = rt_times,y = Intensity,fill="red",alpha=0.6))+theme(legend.position = "none")

  peak_spectra_range <- rt_spectra[rt_box[1]:rt_box[2]]
  peak_rt_range <- rt_times2[rt_box[1]:rt_box[2]]
  # Extract optimal mz range

  # spectrum <- spectra[[which.max(chrom)]]
  spectrum <- spectra[[which.max(chrom3$Intensity)]]
  #plot(spectrum, type = "l")
  section <- (spectrum[spectrum[,1] > (chrom3 %>% filter(Peak_max==T))$m_z-0.3 & spectrum[,1] < (chrom3 %>% filter(Peak_max==T))$m_z+0.3, ])
  # Add a more sophisticated method here
  #plot(section)
  peak_mz_range <- c(min(section[,1]), max(section[,1]))

  subset1 <- map(.x = peaks(mz, peak_spectra_range),.f = function(x) mzfinder2(x, peak_mz_range=peak_mz_range))
  subset2 <- map(.x = peaks(mz, peak_spectra_range),.f = function(x) mzfinder2(x, peak_mz_range=peak_mz_range-1.0078))
  subset3 <- map(.x = peaks(mz, peak_spectra_range),.f = function(x) mzfinder2(x, peak_mz_range=peak_mz_range+1.0078))
  peakmatrix1 <- bind_rows(subset1,.id = "Spectrum")
  peakmatrix2 <- bind_rows(subset2,.id = "Spectrum")
  peakmatrix3 <- bind_rows(subset3,.id = "Spectrum")

  return(chrom3 %>% filter(Peak_max==T) %>% select(D_m_z=m_z,D_rt_times = rt_times,D_rt_times=rt_times) %>%
           mutate(C_mz_low=peak_mz_range[1]
                  ,C_mz_high=peak_mz_range[2]
                  ,C_rt_Start = rt_times2[rt_box[1]]/60
                  ,C_rt_End=rt_times2[rt_box[2]]/60
                  ,B_Area=sum(peakmatrix1$Intensity)
                  ,B_M_1=sum(peakmatrix3$Intensity)
                  ,B_Blank=sum(peakmatrix2$Intensity)) %>%
           select(B_Area,B_Blank,everything()) %>% mutate(
             Peakshape = list(list(chrom3,chrom3c)),C_small_rt_dev=NA)
  )
}


# peakextract_whole plotting ----------------------------------------------
# x <- lipidlist[3,]
# rt=as.numeric(x[1,4])
# mass=as.numeric(x[1,3])
# mz=mz
# ms1_header=ms1_header

# depending on the minimal and maximal time of peaks found:
# REMOVE sometimes, if unused in the future:
plot_whole_dataset_PL <- function(Peaklist, mz, ms1_header,folder,file_name,Date){

  window <- c(500,1000) #lower and upper cut off

  cutoff= 0.1 # cutoff for peakfinder function

  rt_min <-  ms1_header[which.min(abs(ms1_header$retentionTime - min((Peaklist %>% select(lipid,starts_with(match = "C")))$C_rt_Start)*60-12)), 1]
  rt_max <-  ms1_header[which.min(abs(ms1_header$retentionTime - max((Peaklist %>% select(lipid,starts_with(match = "C")))$C_rt_End)*60+12)), 1]

  rt_spectra <- ms1_header$seqNum[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]
  rt_times <- ms1_header$retentionTime[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]

  # Data extraction for the 2D plots
  spectra <- peaks(mz, rt_spectra)
  Spectra <- future_map(spectra,.f = function(x) {tibble(m_z=x[,1],Intensity=x[,2])})
  Spectra <- bind_rows(bind_rows(setNames(Spectra, seq_along(Spectra)), .id = "id")) %>% rename(Spectrum= id)

  # object <- tibble(Spectra=spectra) %>% mutate(Spectra=future_map(Spectra,.f = function(x) {tibble(m_z=x[,1],Intensity=x[,2])}))
  # head(object)
  # chrom_plot <- object %>% unnest(.id = "Spectrum")

  chrom_plot <- left_join(x = Spectra,
                          y = tibble(rt_times2=rt_times/60) %>% rownames_to_column(var = "Spectrum"),
                          by="Spectrum")


  # HIER PASSIERT DAS FILTERING1
  rect_data <- (Peaklist %>% select(lipid,starts_with(match = "C")))
  rect_data_p <- rect_data %>% filter(C_mz_low > window[1] & C_mz_high < window[2])
  plotting <- chrom_plot %>% filter(m_z < max(rect_data_p$C_mz_high,na.rm = T)+1 & m_z > min(rect_data_p$C_mz_low,na.rm = T)-1)
  # plotting$Intensity
  cuter <- 500
  plotting <- plotting %>% filter(Intensity > cuter)
  p <- ggplot()+
    geom_hex(data = plotting, mapping = aes(x = rt_times2,y = m_z,weight=log2(Intensity)),stat = "binhex",binwidth=c(.05,.1))+ # binwidth(x,y)
    scale_fill_gradient2(low = "white",mid = "blue",high = "black",midpoint = max(log10(plotting$Intensity))*10,na.value = "white",)+
    # scale_fill_grey()na.value = "white",low = "yellow",mid = "blue",high = "black",midpoint = 3)+
    geom_rect(data = rect_data_p %>% filter(str_detect(string = lipid,pattern = "PC")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="red",alpha=0.5)+
    geom_rect(data = rect_data_p %>% filter(str_detect(string = lipid,pattern = "PC")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low+1.0078, ymax=C_mz_high+1.0078),fill="red",alpha=0.3)+
    geom_rect(data = rect_data_p %>% filter(str_detect(string = lipid,pattern = "PE")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="blue",alpha=0.5)+
    geom_rect(data = rect_data_p %>% filter(str_detect(string = lipid,pattern = "PE")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low+1.0078, ymax=C_mz_high+1.0078),fill="blue4",alpha=0.3)+
    geom_rect(data = rect_data_p %>% filter(str_detect(string = lipid,pattern = "SM")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="cadetblue1",alpha=0.5)+
    geom_rect(data = rect_data_p %>% filter(str_detect(string = lipid,pattern = "Cer")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="yellow",alpha=0.5)+
    geom_rect(data = rect_data_p %>% filter(str_detect(string = lipid,pattern = "PA")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="chartreuse1",alpha=0.5)+
    geom_rect(data = rect_data_p %>% filter(str_detect(string = lipid,pattern = "PG")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="aquamarine4",alpha=0.5)+
    geom_rect(data = rect_data_p %>% filter(str_detect(string = lipid,pattern = "PI")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="chocolate",alpha=0.5)+
    geom_rect(data = rect_data_p %>% filter(str_detect(string = lipid,pattern = "PS")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="cyan",alpha=0.5)+
    geom_text(data = rect_data_p,mapping = aes(x = (C_rt_End+C_rt_Start)/2,y = C_mz_high,label=lipid),size=2,check_overlap = T)+
    scale_y_continuous(limits = c(min(rect_data_p$C_mz_low,na.rm = T)-1.5,max(rect_data_p$C_mz_high,na.rm = T)+1.5),breaks = seq(from=window[1],to = window[2],by = 2))+
    scale_x_continuous(limits = c(min(rect_data_p$C_rt_Start,na.rm = T)-.1,max(rect_data_p$C_rt_End,na.rm = T)+0.1),breaks = function(limits) seq(from=round(limits[1],digits = 0),to = round(limits[2],digits = 0),by = 0.2))+
    labs(title = unique(Peaklist$Sample))+
    theme_bw()
  # p
  ggsave(filename = here(Folder_created,paste0(file_name,"_PL_",str_c(window,collapse = "-"),"mz_",Sys.Date(),".pdf")),plot = p,device = "pdf",width = 16.5,height = 11.7)
}

# depending on the result4 filtered object only:
plot_whole_dataset_PL2 <- function(result4, mz, ms1_header,file_name,Folder_created){
  Peaklist <- result4 %>% filter(Sample == file_name)
  window <- c(500,1000) #lower and upper cut off

  rt_min <-  ms1_header[which.min(abs(ms1_header$retentionTime - 0*60)), 1]
  rt_max <-  ms1_header[which.min(abs(ms1_header$retentionTime - 12*60+12)), 1]

  rt_spectra <- ms1_header$seqNum[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]
  rt_times <- ms1_header$retentionTime[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]

  # Data extraction for the 2D plots
  spectra <- peaks(mz, rt_spectra)
  Spectra <- future_map(spectra,.f = function(x) {tibble(m_z=x[,1],Intensity=x[,2])})
  Spectra <- bind_rows(bind_rows(setNames(Spectra, seq_along(Spectra)), .id = "id")) %>% rename(Spectrum= id)

  # object <- tibble(Spectra=spectra) %>% mutate(Spectra=future_map(Spectra,.f = function(x) {tibble(m_z=x[,1],Intensity=x[,2])}))
  # head(object)
  # chrom_plot <- object %>% unnest(.id = "Spectrum")

  chrom_plot <- left_join(x = Spectra,
                          y = tibble(rt_times2=rt_times/60) %>% rownames_to_column(var = "Spectrum"),
                          by="Spectrum")

  # HIER PASSIERT DAS FILTERING1
  rect_data <- (Peaklist %>% select(lipid,starts_with(match = "C")))
  rect_data2 <- (Peaklist %>% select(lipid,starts_with(match = "D")))
  rect_data_p <- rect_data %>% filter(C_mz_low > window[1] & C_mz_high < window[2])

  plotting1 <- chrom_plot %>% filter(m_z < 901 & m_z > 499)
  # plotting$Intensity
  cuter <- 100
  plotting <- plotting1 %>% filter(Intensity > cuter)
  p2D <- ggplot()+
    geom_hex(data = plotting, mapping = aes(x = rt_times2,y = m_z,weight=log2(Intensity)),stat = "binhex",binwidth=c(.1,.2))+ # binwidth(x,y)
    # scale_fill_gradient2(low = "white",mid = "blue",high = "black",midpoint = max(log10(plotting$Intensity))*10,na.value = "white")+
    # scale_fill_gradient(low = "white",high = "red",na.value = "white")+
    scale_fill_gradient(low = "white",high = "black",na.value = "white")+
    geom_rect(data = rect_data_p,mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="red",alpha=0.8)+
    geom_point(data = rect_data2,mapping = aes(x=D_rt_times,y = D_m_z),size=.1,alpha=0.8,col="green")+
    geom_text(data = rect_data2,mapping = aes(x=D_rt_times,y = D_m_z+1,label=lipid),size=0.2,alpha=0.8,col="black")+
    # scale_fill_grey()na.value = "white",low = "yellow",mid = "blue",high = "black",midpoint = 3)+
    scale_x_continuous(limits = c(0,12),breaks = function(limits) {seq(from=limits[1],to = limits[2],by = 1)},minor_breaks = function(limits) {seq(from=limits[1],to = limits[2],by = 0.2)},expand = c(0,0))+
    scale_y_continuous(limits = c(500,900),breaks = function(limits) {seq(from=limits[1],to = limits[2],by = 5)},minor_breaks = function(limits) {seq(from=limits[1],to = limits[2],by = 1)},expand = c(0,0))+
    theme_bw()+
    theme(legend.position = "none")

  ggsave(filename = here(Folder_created,paste0(filelist2[i],"_2D_processed_",Sys.Date(),".pdf")),plot = p2D,device = "pdf",width = 16.5,height = 11.7)
}


# as done in the main script for comparability - should write this --------








# tictoc::tic()
# p <- ggplot()+
#   geom_point(data = filtered_data,mapping = aes(x = rt_times2,y = m_z,alpha=Intensity),size=.1)+
#   geom_rect(data = rect_data %>% filter(str_detect(string = lipid,pattern = "PC")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="red")+
#   geom_rect(data = rect_data %>% filter(str_detect(string = lipid,pattern = "PE")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="blue")+
#   geom_rect(data = rect_data %>% filter(str_detect(string = lipid,pattern = "SM")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="brown")+
#   geom_text(data = rect_data,mapping = aes(x = (C_rt_End+C_rt_Start)/2,y = C_mz_high,label=lipid),size=2,check_overlap = T)+
#   labs(title = unique(Peaklist$Sample))+
#   theme_light()
# tictoc::toc()
# ggsave(filename = here(paste0(folder,Date),paste0(file_name,"_PL_A3_",Sys.Date(),".pdf")),plot = p,device = "pdf",width = 16.5,height = 11.7)
# library(plotly)
# ply <- ggplotly(p)
# htmlwidgets::saveWidget(widget = p,file = here(paste0(folder,Sys.Date()),paste0(file_name,"_",Sys.Date(),".html")))




plot_whole_dataset_CL <- function(Peaklist, mz, ms1_header,folder,file_name,Date){

  window <- c(1000,1600) #lower and upper cut off
  cutoff= 0.1 # cutoff for peakfinder function

  rt_min <-  ms1_header[which.min(abs(ms1_header$retentionTime - min((Peaklist %>% select(lipid,starts_with(match = "C")))$C_rt_Start)*60-12)), 1]
  rt_max <-  ms1_header[which.min(abs(ms1_header$retentionTime - max((Peaklist %>% select(lipid,starts_with(match = "C")))$C_rt_End)*60+12)), 1]

  rt_spectra <- ms1_header$seqNum[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]
  rt_times <- ms1_header$retentionTime[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]

  # Data extraction for the 2D plots
  spectra <- peaks(mz, rt_spectra)
  object <- tibble(Spectra=spectra) %>% mutate(Spectra=map(Spectra,.f = function(x) {tibble(m_z=x[,1],Intensity=x[,2])}))

  chrom_plot <- object %>% unnest(.id = "Spectrum")
  chrom_plot <- left_join(x = chrom_plot,
                          y = tibble(rt_times2=rt_times/60) %>% rownames_to_column(var = "Spectrum"),
                          by="Spectrum")

  plotting <- chrom_plot %>% filter(m_z < window[2] & m_z > window[1])
  plotting$Intensity
  cuter <- 100
  plotting <- plotting %>% filter(Intensity > cuter)

  filtered_data <- plotting %>% mutate(rt_times2=round(rt_times2,digits = 2),
                                       m_z=round(x = m_z,digits = 0)) %>%
    group_by(rt_times2,m_z) %>% summarise(Intensity=mean(Intensity))
  rect_data <- (Peaklist %>% select(lipid,starts_with(match = "C")))

  tictoc::tic()
  p <- ggplot()+
    geom_point(data = filtered_data,mapping = aes(x = rt_times2,y = m_z,alpha=Intensity),size=.1)+
    geom_rect(data = rect_data %>% filter(str_detect(string = lipid,pattern = "PC")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="red")+
    geom_rect(data = rect_data %>% filter(str_detect(string = lipid,pattern = "PE")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="blue")+
    geom_rect(data = rect_data %>% filter(str_detect(string = lipid,pattern = "SM")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="brown")+
    geom_text(data = rect_data,mapping = aes(x = (C_rt_End+C_rt_Start)/2,y = C_mz_high,label=lipid),size=2,check_overlap = T)+
    labs(title = unique(Peaklist$Sample))+
    theme_light()
  tictoc::toc()
  p
  ggsave(filename = here(paste0(folder,Date),paste0(file_name,"_CL_A3_",Sys.Date(),".pdf")),plot = p,device = "pdf",width = 16.5,height = 11.7)
  # library(plotly)
  # p <- ggplotly(p)
  # htmlwidgets::saveWidget(widget = p,file = here(paste0(folder,Sys.Date()),paste0(file_name,"_",Sys.Date(),".html")))
}


# generate the point plot, no layer added... ------------------------------

plot_whole_dataset_background <- function(mz, ms1_header){

  window <- c(500,1000) #lower and upper cut off
  cutoff= 0.1 # cutoff for peakfinder function

  rt_min <-  ms1_header[which.min(abs(ms1_header$retentionTime - min((Peaklist %>% select(lipid,starts_with(match = "C")))$C_rt_Start)*60-12)), 1]
  rt_max <-  ms1_header[which.min(abs(ms1_header$retentionTime - max((Peaklist %>% select(lipid,starts_with(match = "C")))$C_rt_End)*60+12)), 1]

  rt_spectra <- ms1_header$seqNum[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]
  rt_times <- ms1_header$retentionTime[ms1_header$seqNum >= rt_min & ms1_header$seqNum <= rt_max]

  # Data extraction for the 2D plots
  spectra <- peaks(mz, rt_spectra)
  object <- tibble(Spectra=spectra) %>% mutate(Spectra=map(Spectra,.f = function(x) {tibble(m_z=x[,1],Intensity=x[,2])}))

  chrom_plot <- object %>% unnest(.id = "Spectrum")
  chrom_plot <- left_join(x = chrom_plot,
                          y = tibble(rt_times2=rt_times/60) %>% rownames_to_column(var = "Spectrum"),
                          by="Spectrum")

  plotting <- chrom_plot %>% filter(m_z < window[2] & m_z > window[1])
  plotting$Intensity
  cuter <- 100
  plotting <- plotting %>% filter(Intensity > cuter)

  filtered_data <- plotting %>% mutate(rt_times2=round(rt_times2,digits = 2),
                                       m_z=round(x = m_z,digits = 0)) %>%
    group_by(rt_times2,m_z) %>% summarise(Intensity=mean(Intensity))
  rect_data <- (Peaklist %>% select(lipid,starts_with(match = "C")))

  tictoc::tic()
  p <- ggplot()+
    geom_point(data = filtered_data,mapping = aes(x = rt_times2,y = m_z,alpha=Intensity),size=1)+
    scale_alpha_continuous(range = c(0,1))+
    geom_rect(data = rect_data %>% filter(str_detect(string = lipid,pattern = "PC")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="red")+
    geom_rect(data = rect_data%>% filter(str_detect(string = lipid,pattern = "PE")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="blue")+
    geom_rect(data = rect_data%>% filter(str_detect(string = lipid,pattern = "SM")),mapping = aes(xmin=C_rt_Start,xmax=C_rt_End,ymin=C_mz_low, ymax=C_mz_high),fill="brown")+
    geom_text(data = rect_data,mapping = aes(x = (C_rt_End+C_rt_Start)/2,y = C_mz_high,label=lipid),size=2,check_overlap = T)+
    labs(title = unique(Peaklist$Sample))
  tictoc::toc()
  p
  ggsave(filename = here(paste0(folder,Sys.Date()),paste0(file_name,"_A3_",Sys.Date(),".pdf")),plot = p,device = "pdf",width = 16.5,height = 11.7)
  library(plotly)
  p <- ggplotly(p)
  htmlwidgets::saveWidget(widget = p,file = here(paste0(folder,Sys.Date()),paste0(file_name,"_",Sys.Date(),".html")))
}



# peaklist_regression_cal_correction plot ---------------------------------
ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}