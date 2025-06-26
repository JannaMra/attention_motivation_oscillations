# EDA Analyse
# Design:
# cue 1: high intensity shock (shock)
# cue 2: high intensity flight (flight)
# cue 3: low intensity shock (safety)
# cue 4: low intensity flight (active control)

#packages 
#library(signal)
library(tidyverse)
library(stringr)
library(apa)
library(ggbeeswarm)
library(afex)
options(dplyr.summarise.inform = FALSE)

path <- paste0(getwd(),"/")


partial_eta_squared_ci <- function(x)
{
  # Calculate upper and lower limit of Lambda for 90%-CI
  limits <- MBESS::conf.limits.ncf(x$F, .95, x$DFn, x$DFd)
  
  # Convert Lambda to partial eta-squared
  ci <- c(limits$Lower.Limit, limits$Upper.Limit) /
    (c(limits$Lower.Limit, limits$Upper.Limit) + x$DFn + x$DFd + 1)
  
  # Replace NA with zero (`conf.limits.ncf` returns NA instead of 0)
  ci[is.na(ci)] <- 0
  
  ci <- sprintf("%.2f", round(ci, 2))
  
  paste0("[", ci[1], "; ", ci[2], "]")
}

#global parameters, options and plots

ucr_plots = T
cr_plots = T

downsampling = F
sample_rate_new = 20 #for downsampling
check_ds_plots = F

lowpass = F
lowPassFreq = 10 #low pass filter (Hz)
highpass = F
highPassFreq = 0.01 #high pass filter (Hz)

check_minmax_plots = F 
saveSclPlotsToFile = T


shock_triggers = c(5,6,7,8) #Vector with the names of the US triggers

sample_rate = 20  #samplerate after export
trial_length = 12  #length of trial in s

#Min-Max-SCR Analysis:

ucrMinWindow = c(.5, 3) #within how many seconds after UCS must the detected minimum appear?
ucrMin = .1  #minimum amplitude of an individual unconditioned response (in µS)
ucrRiseMax = 6 #latest Peak of UCR in s

crMin = 0.02      #minimum amplitude of CS-SCR
crMinWindow = c(.8,3) #within how many seconds after CS must the detected minimum appear? (Bouton et al., 2012)
crRiseMax = 6 #latest Peak of CR in s (Bouton et al., 2012)

#SCL Analysis:

scl_bins = T #calculate SCL in 1s bins after CS onset
baselineWindow = c(-.5,0) #correct for Baseline in this time window
sclWindow = c(0,10) #calculate SCL in this thime window


#First and second responses (FIR & SIR) analysis: 

fir_algorithm = F

eirWindow = c(1,10) #entire interval response window
firWindow = c(1,5) #first interval response window
sirWindow = c(5,10) #second interval response window
eirMinWindow = c(-.5,3) #minimun interval window for baseline correction

{ # Functions ---------------------------------------------------------------
  sample.down = function(signal, conversionRate) {
    suppressWarnings(matrix(signal, ncol=conversionRate, byrow=T)) %>% apply(1, mean) %>% 
      head(-1) #discard last sample because it gets distorted by zero padding
  }
  
  closestRepresentative = function(x, values, returnIndices=F) {
    indices = x %>% sapply(function(x, data) which.min(abs(data-x)), data=values)
    if (returnIndices) return(indices)
    return(values[indices])
  }
  
  localMaxima = function(signal, integerMode=F) {
    result = signal %>% c(-ifelse(integerMode, .Machine$integer.max, Inf), .) %>% #prepend smallest value to allow that first sample point is also local maximum
      diff() %>% {. > 0} %>% #which sample points are rising compared to their predecessor?
      rle() %>% .$lengths %>% #for how long are segments rising consecutively?
      cumsum() %>% #cumulative sum to find indices of local extrema
      {.[seq.int(from=1, to=length(.), by=2)]} #only take every other point
    if (signal[[1]] == signal[[2]]) result = result[-1]
    
    return(result)
  }
  
  localMinima = function(signal, integerMode=F) {
    result = signal %>% c(ifelse(integerMode, .Machine$integer.max, Inf), .) %>% #prepend biggest value to allow that first sample point is also local maximum
      diff() %>% {. > 0} %>% #which sample points are rising compared to their predecessor?
      rle() %>% .$lengths %>% #for how long are segments rising consecutively?
      cumsum() %>% #cumulative sum to find indices of local extrema
      {.[seq.int(from=1, to=length(.), by=2)]} #only take every other point
    if (signal[[1]] == signal[[2]]) result = result[-1]
    
    return(result)
  }
  
  inflectionPoints = function(signal, integerMode=F, upOnly=T) {
    extrema = signal %>% c(ifelse(integerMode, .Machine$integer.max, Inf), .) %>% diff()
    curvature = extrema %>% c(ifelse(integerMode, .Machine$integer.max, Inf)) %>% #for curvature, append biggest value (first and last sample point have no defined curvature)
      diff()
    result = curvature %>% #{. == 0} %>% which() #not working for discrete graphs
      {. < 0} %>% #which sample points are curved downwards (compared to adjacent)?
      rle() #run length of consecutive curvature
    
    if (upOnly) keep = result$values
    
    result = result %>% .$lengths %>% #for how long are segments curved upwards consecutively?
      cumsum() #cumulative sum to find indices of inflection points
    
    if (upOnly) result = result[keep]
    
    return(result)
  }
}


edas_cr_list = list()
edas_list = list()
edas_df_list = list()
eda_maxima_list = list()
eda_minima_list = list()
eda_inflection_list = list()
eda_df = tibble()
eda_unified = tibble()

filemat = list.files("C:/Users/jat41gk/Documents/Projekte/Action-Motivation-Oscillations/Analyse/Data/Physio/EXP/", pattern=".*_EDA.txt") # Das sollte der Ordner sein, in dem deine ganzen Files liegen. Pro VP ein File.
filemat <- filemat[!grepl("vp05|vp22|vp30|vp51", filemat)]  #filter general exclusions


eda_ucr = tibble(ID = filemat %>% sub("_EDA.txt","", .), ucr = 0, valid = 0, lat = 0, rise = 0) %>% 
  #separate(subject,c("Exp","ID","tech","tech2"),sep="_",remove=T) %>%
  select(ID,ucr,valid,lat,rise)
eda_cr = tibble(ID = filemat %>% sub("_EDA.txt","", .), cr = 0, valid = 0, lat = 0, rise = 0) %>% 
  #separate(subject,c("Exp","ID","tech","tech2"),sep="_",remove=T) %>%
  select(ID,cr,valid,lat,rise)



for (subject_inmat in filemat){ #Jetzt rechnet er den Spaß für jedes File durch, wenn du einzelne Files berechnen willst, mach es nicht mit der for loop sondern kommentier die nächste zeile wieder ein
  
  #subject_inmat = filemat[[1]]
  eda <- read.table(paste0("C:/Users/jat41gk/Documents/Projekte/Action-Motivation-Oscillations/Analyse/Data/Physio/EXP/",subject_inmat), sep="\t",stringsAsFactors =F) %>% #read export
    mutate(V1= str_replace(gsub("\\s+", " ", str_trim(V1)), "B", "b")) %>% separate(V1,c("time","EDA","trigger"),sep=" ") %>%   #remove whitespace and extraxt 3 channels
    mutate(sample = 1:n(),time=suppressWarnings(as.numeric(time)), EDA=suppressWarnings(as.numeric(EDA)),trigger=suppressWarnings(as.numeric(trigger))) %>% #channels as numeric
    mutate(EDA = -EDA)
  
  filename =  filemat[filemat == subject_inmat] %>% str_remove(.,pattern="_EDA.txt")# %>% str_remove(.,pattern="SiCP_")
  
  for (row in 1:(nrow(eda)-1)){ if (eda$trigger[row] + eda$trigger[row+1] == 14) {eda$trigger[row+1]=0}} #remove extra shock-markers
  
  #sum(eda$trigger == 2);sum(eda$trigger == 3);sum(eda$trigger == 4);sum(eda$trigger == 5);sum(eda$trigger == 6);sum(eda$trigger == 7);
  #sum(eda$trigger == 8);sum(eda$trigger == 9);sum(eda$trigger == 10);sum(eda$trigger == 11);sum(eda$trigger == 12); #count triggers and check
  
  
  # if (filename == "tfo06"){ #VP06 restarted the experiment after 16 trials because of technical issues.
  #   eda <- eda %>% filter(time > 440)  #Remove first 440s from data
  #   
  # }

  if (sum(eda$trigger == 1) != 40) { 
    warning(paste0("Warning: Too few/many Threat triggers"))
    eda$trigger[head(which(eda$trigger == 1), sum(eda$trigger == 1) -40)] <- 0}
    if (sum(eda$trigger  == 2) != 40) { 
      warning(paste0("Warning: Too few/many  Flight triggers"))
      eda$trigger[head(which(eda$trigger == 2), sum(eda$trigger == 2) -40)] <- 0}
    if (sum(eda$trigger  == 3) != 40) { 
      warning(paste0("Warning: Too few/many  Safe triggers"))
      eda$trigger[head(which(eda$trigger == 3), sum(eda$trigger == 3) -40)] <- 0}
    if (sum(eda$trigger  == 4) != 40) { 
      warning(paste0("Warning: Too few/many  Safe triggers"))
      eda$trigger[head(which(eda$trigger == 4), sum(eda$trigger == 4) -40)] <- 0}
  
  
  
  print("Data reading complete!")
  
  
  if (downsampling) {
    #downsample (all columns)
    triggers_time = eda$time[eda$trigger != 0] #eda$time[eda$Trigger %>% is.na() == FALSE]
    conversion = round(sample_rate / sample_rate_new)
    eda_downsampled = data.frame(time=sample.down(eda$time,conversion),
                                 EDA=sample.down(eda$EDA,conversion),
                                 trigger=0) %>% 
      mutate(sample=1:n()) %>% select(sample, everything())
    
    #calculate closest position for trigger onsets in downsampled time
    triggers_time_old = eda$time[eda$trigger != 0]
    triggers_indices_new = triggers_time_old %>% closestRepresentative(eda_downsampled$time, returnIndices = T) #for each old trigger time, find index of closest existing downsampled time
    eda_downsampled$trigger[triggers_indices_new] = eda$trigger[eda$trigger != 0] #inject conditions as triggers of onsets
    
    # eda_downsampled %>% ggplot(aes(x=time, y=eda)) + geom_line() + geom_vline(xintercept = eda$time[eda$trigger!=0])
    if(check_ds_plots){
      eda %>% ggplot(aes(x=time, y=EDA)) + geom_line() + geom_line(data=eda_downsampled, color="red",linetype=4) + 
        geom_vline(xintercept = eda$time[eda$trigger!=0],alpha=0.1)  %>% print()
    }  
    
    eda = eda_downsampled; rm(eda_downsampled)
    print("Downsampling complete!")
    
  }  
  
  #smoothing / filtering
  if (lowpass) {
    reps = 100 #filters can produce artifacts on the edges => fill edges with first and last values
    eda_filtered = c(rep(first(eda$EDA), reps), eda$EDA, rep(last(eda$EDA), reps)) %>% 
      signal::filtfilt(signal::butter(2, lowPassFreq/(ifelse(downsampling, sample_rate_new, sample_rate)/2)), .) #low-pass filter
    
    if(highpass){
      eda_filtered = eda_filtered %>%  # high-pass filter
        signal::filtfilt(signal::butter(2, type="high", highPassFreq/(ifelse(downsampling, sample_rate_new, sample_rate)/2)), .)  
    }
    
    eda_filt = eda
    eda_filt$EDA = eda_filtered[(reps+1):(length(eda_filtered)-reps)]; rm(eda_filtered)
    
    if (check_smoothed_plots) {
      eda %>% ggplot(aes(x=time, y=EDA)) + geom_line() + geom_line(data=eda_filt, color="red",linetype=4) + 
        geom_vline(xintercept = eda$time[eda$trigger!=0],alpha=0.1)  %>% print()
    }
    print("lowpass filtering complete!")
    
  }  
  
  
  eda = eda %>% mutate(us = case_when(trigger %in% shock_triggers ~ T,
                                      TRUE ~ F))
  eda_vp = eda %>% dplyr::filter(trigger != 0) %>% select(-EDA) %>% #One trial per row per VP
    mutate(trial = 1:n(), condition = trigger,
           us = us, us_prior = c(FALSE, lag(us)[-1]),
           time_start = time, time_end = time + trial_length,
           sample_start = sample, 
           sample_end =  {sample+(trial_length*ifelse(downsampling,sample_rate_new,sample_rate))} %>% round()) %>% 
    select(-trigger, -time, -sample) %>% select(trial, condition, everything())
  
  
  edas_df_list[[filename]] = eda_vp
  edas_list[[filename]] = eda
  
  # no nesting anymore
  #eda_vp$eda = eda_vp$time_start %>% lapply(function(start) 
  #  eda %>% select(-us) %>%  filter(time >= start, time < start + trial_length)) 
  
  eda_maxima = tibble(sample = localMaxima(eda$EDA)) %>% 
    mutate(n = 1:n(), time=eda$time[sample], EDA= eda$EDA[sample]) %>% select(n, everything())
  eda_minima = tibble(sample = localMinima(eda$EDA)) %>% 
    mutate(n = 1:n(), time=eda$time[sample], EDA = eda$EDA[sample]) %>% select(n, everything())
  eda_inflection = tibble(sample = inflectionPoints(eda$EDA)) %>% 
    mutate(n = 1:n(), time=eda$time[sample], EDA = eda$EDA[sample]) %>% select(n, everything())
  
  eda_maxima_list[[filename]] = eda_maxima
  eda_minima_list[[filename]] = eda_minima
  eda_inflection_list[[filename]] = eda_inflection
  
  if (check_minmax_plots || saveSclPlotsToFile) {
    sclPlot = eda %>% ggplot(aes(x=time, y=EDA)) + 
      geom_vline(xintercept = eda$time[eda$trigger!=0], color=ifelse(eda_vp$us, "orange", "black")) +
      geom_line() + 
      #geom_point(data=eda.inflection, color="purple") + 
      geom_point(data=eda_maxima, color="red") + 
      geom_point(data=eda_minima, color="blue") + 
      ggtitle(filename)
    if (check_minmax_plots) sclPlot %>% print()
    if (saveSclPlotsToFile) sclPlot %>% ggsave(paste0("../Plots/EDA/", filename, ".png"), plot=., device="png", width=1920/50, height=1080/300, dpi=300)
  }
  
  # UCR Scoring    
  eda_vp_ucr = eda_vp %>% filter(us == TRUE)
  
  for (t in 1:nrow(eda_vp_ucr)) {
    
    #t = 1
    start = eda_vp_ucr$time_start[t]
    eda_minima_trial = eda_minima %>% filter(time >= start + min(ucrMinWindow),
                                             time <= start + max(ucrMinWindow))
    eda_maxima_trial = eda_maxima %>% filter(time > suppressWarnings(min(eda_minima_trial$time))) %>% 
      head(nrow(eda_minima_trial)) #same amount of maxima as minima
    
    #if no minimum found within time range use deflection point as minimum
    if (nrow(eda_minima_trial)==0) { 
      eda_minima_trial = eda_inflection %>% filter(time >= start + min(ucrMinWindow),
                                                   time <= start + max(ucrMinWindow))
      eda_maxima_trial = eda_maxima %>% filter(time > suppressWarnings(min(eda_minima_trial$time))) %>% 
        head(nrow(eda_minima_trial)) #same amount of maxima as minima
      
      while (nrow(eda_maxima_trial)>1) {
        latestMaxTimes = eda_maxima_trial %>% tail(2) %>% .$time
        if (eda_minima_trial %>% filter(time < max(latestMaxTimes), #check if any inflection point is between last two maxima (may be false if two inflection points before earlier max)
                                        time > min(latestMaxTimes)) %>% 
            nrow() %>% {. > 0}) break
        else { #if no inflection point between last two maxima
          eda_minima_trial = eda_minima_trial %>% head(-1) #delete last inflection point
          eda_maxima_trial = eda_maxima_trial %>% head(nrow(eda_minima_trial)) #same amount of maxima as minima
        }
      }
      
      #check if true min occurs between inflection point and corresponding max (with true min being too late)
      toKeep = c()
      if (nrow(eda_minima_trial) > 0) { 
        for (i in 1:nrow(eda_minima_trial)) {
          if (eda_minima %>% filter(time > eda_minima_trial[i, "time"], time < eda_maxima_trial[i, "time"]) %>% 
              nrow() %>% {. == 0}) toKeep = toKeep %>% c(i)
        }
        eda_minima_trial = eda_minima_trial[toKeep,]; eda_maxima_trial = eda_maxima_trial[toKeep,] #only keep inflection points that don't have a true minimum between themself and their max
      }
    }
    
    #calculate SCRs
    eda_minima_trial_renamed = eda_minima_trial; names(eda_minima_trial_renamed) = names(eda_minima_trial_renamed) %>% paste0(".min")
    eda_maxima_trial_renamed = eda_maxima_trial; names(eda_maxima_trial_renamed) = names(eda_maxima_trial_renamed) %>% paste0(".max")
    scr = bind_cols(eda_minima_trial_renamed, eda_maxima_trial_renamed)
    #names(scr) = names(scr) %>% {ifelse(endsWith(., "1"), gsub("1", ".max", ., fixed=T), paste(., "min", sep="."))}
    scr = scr %>% mutate(SCR = EDA.max - EDA.min,
                         Lat = time.min - start,
                         Rise = time.max - start) %>% select(SCR, Lat, Rise, everything()) %>% 
      filter(SCR == max(SCR)) #take maximum SCR
    
    if (nrow(scr) > 1){scr <- scr[1,]}
    
    if (nrow(scr) == 0) {scr = scr %>% bind_rows(tibble(SCR = 0, Lat = NA, Rise = NA)) #if no valid SCR -> 0
    } else if (scr$SCR <= ucrMin) {scr = scr %>% mutate(SCR = 0, Lat = NA, Rise = NA)
    } else if (scr$Rise > ucrRiseMax ) {scr = scr %>% mutate(SCR = 0, Lat = NA, Rise = NA)}
    
    scr_out = scr %>% select(-c(n.min,sample.min,n.max,sample.max))
    eda_vp_ucr[t, names(scr_out)] = scr_out
    
  }  
  
  eda_ucr[eda_ucr$ID==filename, -1] = eda_vp_ucr %>% mutate(valid = SCR > ucrMin) %>% 
    summarise(scr = mean(SCR, na.rm=T), 
              valid = mean(valid),
              lat = mean(Lat, na.rm=T),
              rise = mean(Rise, na.rm=T))
  
  
  if (ucr_plots) {
    
    unified = eda_vp_ucr %>% mutate(EDA = time_start %>% lapply(function(start) 
      eda %>% filter(time >= start,
                     time <= start + 8) %>% select(-trigger))
    )
    
    for (t in 1:nrow(unified)) {
      unified$EDA[[t]] = unified$EDA[[t]] %>% 
        mutate(trial = t, sample = sample - min(sample) + 1, time = time - min(time)) #unify starting time to allow overlap
      
      
    }
    
    minmax_unified = unified %>% select(trial,EDA.min,time.min,EDA.max,time.max,time_start) %>% 
      mutate(time.min = time.min - time_start,time.max = time.max - time_start,trial = as.factor(1:n())) %>% 
      filter(time.max < 8)
    
    
    unified$EDA %>% bind_rows() %>% 
      mutate(trial = as.factor(trial)) %>% 
      {ggplot(., aes(x=time, y=EDA, color=trial,fill=trial)) + 
          geom_vline(xintercept=ucrMinWindow, color="grey") + #borders of min/max scoring
          geom_path() + scale_color_viridis_d() + scale_fill_viridis_d() +
          geom_point(data=minmax_unified,aes(x=time.min,y=EDA.min),shape=25,size=1)+
          geom_point(data=minmax_unified,aes(x=time.max,y=EDA.max),shape=24,size=1)+
          ggtitle(filename) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))} %>% 
      #print()
      ggsave(paste0("C:/Users/jat41gk/Documents/Projekte/Action-Motivation-Oscillations/Analyse/Data/Physio/EDA_plots/ucr/", filename, ".png"), plot=., device="png", width=1920/150, height=1080/150, dpi=300)
    
    
    
  }
  
  print(paste0("UCRs: ", filename, " of ", length(filemat), " files processed!"))
  
  eda_vp_cr = eda_vp %>% filter(us == F)
  
  for (t in 1:nrow(eda_vp_cr)) {
    
    #t <- 2
    start = eda_vp_cr$time_start[t]
    
    baseline_trial = eda %>% filter(time <= start + max(baselineWindow),
                                    time >= start + min(baselineWindow)) %>% 
      summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
    
    
    scl_trial = eda %>% filter(time <= start + max(sclWindow),
                               time >= start + min(sclWindow)) %>% 
      summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
 
    if(scl_bins){
      
      scl_1 <- eda %>% filter(time >= start,
                              time <= start+1) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_2 <- eda %>% filter(time >= start+1,
                              time <= start+2) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_3 <- eda %>% filter(time >= start+2,
                              time <= start+3) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_4 <- eda %>% filter(time >= start+3,
                              time <= start+4) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_5 <- eda %>% filter(time >= start+4,
                              time <= start+5) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_6 <- eda %>% filter(time >= start+5,
                              time <= start+6) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_7 <- eda %>% filter(time >= start+6,
                              time <= start+7) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_8 <- eda %>% filter(time >= start+7,
                              time <= start+8) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_9 <- eda %>% filter(time >= start+8,
                              time <= start+9) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_10 <- eda %>% filter(time >= start+9,
                              time <= start+10) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_11 <- eda %>% filter(time >= start+10,
                               time <= start+11) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_12 <- eda %>% filter(time >= start+11,
                               time <= start+12) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_13 <- eda %>% filter(time >= start+12,
                               time <= start+13) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      scl_14 <- eda %>% filter(time >= start+13,
                               time <= start+14) %>%
        summarize(mean = mean(EDA,na.rm=T)) %>% unlist()
      
      }
     
    if(fir_algorithm){
    
    eir_max_searching = T
    upperLimit <- start + max(eirWindow)
    lowerLimit <- start + min(eirWindow)
    
    while (eir_max_searching){ 
      
      eir_max = eda %>% 
        filter(time <= upperLimit &
                 time >= lowerLimit) %>%
        filter(EDA == max(EDA)) %>% head(1)
      
      
      if(eir_max$time > lowerLimit & eir_max$time < upperLimit){
        eir_max_searching <- F
      }    
      
      if (isTRUE(all.equal(eir_max$time,upperLimit))){
        upperLimit <- round(upperLimit - 1/sample_rate,2)
      } else if(isTRUE(all.equal(eir_max$time,lowerLimit))){
        lowerLimit <- round(lowerLimit + 1/sample_rate,2)
      }
      
      if(upperLimit == lowerLimit){
        eir_max$EDA <- NA
        eir_max_searching <- F
      }
    }
    
    
    
    fir_max_searching = T
    upperLimit <- start + max(firWindow)
    lowerLimit <- start + min(firWindow)
    
    while (fir_max_searching){ 
      
      fir_max = eda %>% 
        filter(time <= upperLimit &
                 time >= lowerLimit) %>%
        filter(EDA == max(EDA)) %>% head(1)
      
      
      if(fir_max$time > lowerLimit & fir_max$time < upperLimit){
        fir_max_searching <- F
      }    
      
      if (isTRUE(all.equal(fir_max$time,upperLimit))){
        upperLimit <- round(upperLimit - 1/sample_rate,2)
      } else if(isTRUE(all.equal(fir_max$time,lowerLimit))){
        lowerLimit <- round(lowerLimit + 1/sample_rate,2)
      }
      
      if(upperLimit == lowerLimit){
        fir_max$EDA <- NA
        fir_max_searching <- F
      }
    }
    
    
    sir_max_searching = T
    upperLimit <- start + max(sirWindow)
    lowerLimit <- start + min(sirWindow)
    
    while (sir_max_searching){ 
      
      sir_max = eda %>% 
        filter(time <= upperLimit &
                 time >= lowerLimit) %>%
        filter(EDA == max(EDA)) %>% head(1)
      
      
      if(sir_max$time > lowerLimit & sir_max$time < upperLimit){
        sir_max_searching <- F
      }    
      
      if (isTRUE(all.equal(sir_max$time,upperLimit))){
        upperLimit <- round(upperLimit - 1/sample_rate,2)
      } else if(isTRUE(all.equal(sir_max$time,lowerLimit))){
        lowerLimit <- round(lowerLimit + 1/sample_rate,2)
      }
      
      if(upperLimit == lowerLimit){
        sir_max$EDA <- NA
        sir_max_searching <- F
      }
    }
    
    eir_min = eda %>% 
      filter(time <= start + max(eirMinWindow) &
               time >= start + min(eirMinWindow)) %>%
      filter(EDA == min(EDA)) %>% head(1)
    
    }
    
    eda_minima_trial = eda_minima %>% filter(time >= start + min(crMinWindow),
                                             time <= start + max(crMinWindow)) 
    eda_maxima_trial = eda_maxima %>% filter(time > suppressWarnings(min(eda_minima_trial$time))) %>% 
      head(nrow(eda_minima_trial)) #same amount of maxima as minima
    
    if (nrow(eda_maxima_trial) < nrow(eda_minima_trial)) { #if recording ended before next max
      eda_maxima_trial = eda_maxima_trial %>% bind_rows(
        eda %>% filter(time >= max(eda_minima_trial$time)) %>% filter(EDA==max(EDA)) %>% mutate(n = 0) %>% select(n, everything()) %>% select(-trigger)
      )
    }
    
    
    #if no minimum found within time range use deflection point as minimum
    if (nrow(eda_minima_trial)==0) { 
      eda_minima_trial = eda_inflection %>% filter(time >= start + min(crMinWindow),
                                                   time <= start + max(crMinWindow))
      eda_maxima_trial = eda_maxima %>% filter(time > suppressWarnings(min(eda_minima_trial$time))) %>% 
        head(nrow(eda_minima_trial)) #same amount of maxima as minima
      if (nrow(eda_maxima_trial) < nrow(eda_minima_trial)) { #if recording ended before next max
        eda_maxima_trial = eda_maxima_trial %>% bind_rows(
          eda %>% filter(time >= max(eda_minima_trial$time)) %>% filter(EDA==max(EDA)) %>% mutate(n = 0) %>% select(n, everything()) %>% select(-trigger)
        )
      }
      
      #check if mins and maxs alternate (in case inflection points were used)
      while (nrow(eda_maxima_trial)>1) {
        latestMaxTimes = eda_maxima_trial %>% tail(2) %>% .$time
        if (eda_minima_trial %>% filter(time < max(latestMaxTimes), #check if any inflection point is between last two maxima (may be false if two inflection points before earlier max)
                                        time > min(latestMaxTimes)) %>% 
            nrow() %>% {. > 0}) break
        else { #if no inflection point between last two maxima
          eda_minima_trial = eda_minima_trial %>% head(-1) #delete last inflection point
          eda_maxima_trial = eda_maxima_trial %>% head(nrow(eda_minima_trial)) #same amount of maxima as minima
        }
      }
      
      #check if true min occurs between inflection point and corresponding max (with true min being too late)
      toKeep = c()
      if (nrow(eda_minima_trial) > 0) { 
        for (i in 1:nrow(eda_minima_trial)) {
          if (eda_minima %>% filter(time > eda_minima_trial[i, "time"], time < eda_maxima_trial[i, "time"]) %>% 
              nrow() %>% {. == 0}) toKeep = toKeep %>% c(i)
        }
        eda_minima_trial = eda_minima_trial[toKeep,]; eda_maxima_trial = eda_maxima_trial[toKeep,] #only keep inflection points that don't have a true minimum between themself and their max
      }
    }
    
    #calculate SCRs
    eda_minima_trial_renamed = eda_minima_trial; names(eda_minima_trial_renamed) = names(eda_minima_trial_renamed) %>% paste0(".min")
    eda_maxima_trial_renamed = eda_maxima_trial; names(eda_maxima_trial_renamed) = names(eda_maxima_trial_renamed) %>% paste0(".max")
    scr = bind_cols(eda_minima_trial_renamed, eda_maxima_trial_renamed)
    #names(scr) = names(scr) %>% {ifelse(endsWith(., "1"), gsub("1", ".max", ., fixed=T), paste(., "min", sep="."))}
    scr = scr %>% mutate(SCR = EDA.max - EDA.min,
                         Lat = time.min - start,
                         Rise = time.max - start,
                         Baseline = baseline_trial,
                         SCL = scl_trial) %>% select(SCR, Lat, Rise, SCL, Baseline, everything()) %>% 
      filter(SCR == max(SCR)) %>% #take maximum SCR
      head(1) #if two or more equal SCRs, take only the first
    
    if (nrow(scr) == 0) {scr = scr %>% bind_rows(tibble(SCR = 0, Lat = NA, Rise = NA, Baseline = baseline_trial, SCL = scl_trial))
    } else if (scr$SCR <= crMin) {scr = scr %>% mutate(SCR = 0, Lat = NA, Rise = NA, Baseline = baseline_trial, SCL = scl_trial)
    } else if (scr$Rise > crRiseMax) {scr = scr %>% mutate(SCR = 0, Lat = NA, Rise = NA, Baseline = baseline_trial, SCL = scl_trial)}
    
    
    scr_out = scr %>% select(1:5,8,9,12,13)

    if(scl_bins){
      scr_out <- scr_out %>%
        mutate(
          scl_1 = scl_1 - Baseline,
          scl_2 = scl_2 - Baseline,
          scl_3 = scl_3 - Baseline,
          scl_4 = scl_4 - Baseline,
          scl_5 = scl_5 - Baseline,
          scl_6 = scl_6 - Baseline,
          scl_7 = scl_7 - Baseline,
          scl_8 = scl_8 - Baseline,
          scl_9 = scl_9 - Baseline,
          scl_10 = scl_10 - Baseline,
          scl_11 = scl_11 - Baseline,
          scl_12 = scl_12 - Baseline,
          scl_13 = scl_13 - Baseline,
          scl_14 = scl_14 - Baseline)
    }
    
    
    
    if(fir_algorithm){
      scr_out = scr %>% select(1:5,8,9,12,13) %>%
        mutate(EDA.eir = eir_max$EDA,
               EDA.eir.lat = eir_max$time - start,
               EDA.fir = fir_max$EDA,
               EDA.fir.lat = fir_max$time - start,
               EDA.sir = sir_max$EDA,
               EDA.sir.lat = sir_max$time - start,
               EDA.eir.min = eir_min$EDA,
               EDA.eir.min.lat = eir_min$time - start
        )
    }
    
    eda_vp_cr[t, names(scr_out)] = scr_out
    #readline(prompt="Press [enter] to continue")
  } 
  
  eda_vp_cr = eda_vp_cr %>% group_by(condition) %>% mutate(trial_condition = row_number()) %>% ungroup()%>% mutate(trial = 1:160)
  edas_cr_list[[filename]] = eda_vp_cr
  eda_df  = rbind(eda_df,eda_vp_cr %>% mutate(ln_cr = log(SCR+1),ID=filename)) 
 
  print("trialwise parameterization complete!")
  
  if (cr_plots) {
    
    unified = eda_vp_cr %>% mutate(EDA = time_start %>% lapply(function(start) 
      eda %>% filter(time >= start,
                     time <= start + trial_length) %>% select(-trigger)))
    
    
    for (t in 1:nrow(unified)) {
      unified$EDA[[t]] = unified$EDA[[t]] %>% 
        mutate(trial = t, time = time - min(time), condition = unified$condition[[t]],
               trial_condition = unified$trial_condition[[t]]) #unify starting time to allow overlap
    }
    
    
    minmax_unified = unified %>% select(trial,condition,EDA.min,time.min,EDA.max,time.max,time_start) %>% 
      mutate(time.min = time.min - time_start,time.max = time.max - time_start,trial = as.factor(1:n())) %>% 
      filter(time.max < 6)
    
    
    unified$EDA %>% bind_rows() %>% 
      mutate(trial = as.factor(trial)) %>% 
      {ggplot(., aes(x=time, y=EDA, color=trial,fill=trial)) + facet_wrap(vars(condition)) +
          geom_vline(xintercept=c(ucrMinWindow,crRiseMax), color="blue",linetype="dashed") + #borders of min/max scoring
          geom_path() + scale_color_viridis_d() + scale_fill_viridis_d() +
          geom_point(data=minmax_unified,aes(x=time.min,y=EDA.min),shape=25,size=1)+ 
          geom_point(data=minmax_unified,aes(x=time.max,y=EDA.max),shape=24,size=1)+ 
          ggtitle(filename) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))} %>% 
      #print()
      ggsave(paste0("C:/Users/jat41gk/Documents/Projekte/Action-Motivation-Oscillations/Analyse/Data/Physio/EDA_plots/cr/", filename, ".png"), plot=., 
             device="png", width=1920/150, height=1080/150, dpi=300)
    
    print("CR plotting complete!")
  }
  
  unified <- unified %>% mutate(ID = filename) %>% mutate(trial = 1:160)
  eda_unified = rbind(eda_unified,unified)
  
  
  print(paste0("CRs: ", filename, " of ", length(filemat), " files processed!"))
  
}#end inmat for loop 

for (t in 1:nrow(eda_unified)) {
  eda_unified$EDA[[t]] = eda_unified$EDA[[t]] %>% 
    mutate(EDA_bl = EDA - eda_unified$Baseline[[t]], ID = eda_unified$ID[[t]]) #unify starting time to allow overlap
}

#rm(eda, eda_vp)

# Read or Save ga_unified and pupil_df for further processing ------------------

# saveRDS(eda_unified,"EDA_unified.RData")
# saveRDS(eda_df,"EDA_df.RData")

#######################################################
## FOR DATA ANALYSIS READ PREPROCESSED FILES HERE #####
#######################################################

eda_unified <- readRDS("EDA_unified.RData")
eda_df <- readRDS("EDA_df.Rdata")

# Inference statistics-------------------------------------------------------------

# read in prot files 
vpn <- paste("vp",ifelse(1:52<10,"0",""),1:52,sep="")
vpn <- vpn[!(vpn %in% exclusions)]  # Exclusions

All_Prot_Data <- data.frame()

for (vp in vpn){
  
  print(vp)
  
  All_Prot_Data_VP <- read.csv2(paste(path.prot,vp,".csv",sep=""), header = TRUE, sep=";", dec=",", check.names=F)
  All_Prot_Data_VP$vp <- vp
  
  All_Prot_Data <- rbind(All_Prot_Data, All_Prot_Data_VP)
  
}

write.csv2(All_Prot_Data,paste0(path.prot,"All_Prots.csv"),row.names=FALSE)
All_Prot_Data <- read.csv2(paste0(path.prot,"All_Prots.csv"))

# filter out problematic trials

eda <- eda_df %>% rename(vp = ID)
eda <- full_join(All_Prot_Data, eda)

eda <- eda %>%
  filter(problem == 0)


eda_long <- eda %>%  #filter(!(ID %in% c("tfo07"))) %>% 
  select(vp,scl_1:scl_14,condition,trial, trial_condition) %>%
  pivot_longer(scl_1:scl_14, names_to = "timebin", values_to ="scl") %>%
  separate(timebin,c("quark","timebin")) %>% mutate(timebin = as.numeric(timebin))  %>% 
  #filter(trial_condition <= 20) %>%
  select(-quark) 

eda_long <- eda_long %>% 
  filter(timebin <= 10) %>% 
  #filter(timebin >= 3 & timebin <= 10) %>% #filter(valid == T) %>%
  group_by(vp,condition, timebin) %>%  summarise(scl = mean(scl)) %>%
  mutate(vp = as.factor(vp), condition = as.factor(condition),
         ActPass = case_when(condition ==1 | condition == 3 ~ "Passive",
                             condition == 2 | condition == 4 ~ "Active"),
         Intensity = case_when(condition == 1 | condition == 2 ~ "High",
                             condition== 3 | condition == 4 ~ "Low")) 

# saveRDS(eda_long,"EDA_long.RData")

eda_long %>% 
  ez::ezANOVA(dv=.(scl), wid=.(vp), within=.c(ActPass,Intensity), 
              detailed=T, type=2) %>% anova_apa()


anova <- eda_long %>% 
  group_by(vp,ActPass, Intensity,timebin) %>%  summarise(scl = mean(scl)) %>%
  mutate(vp = as.factor(vp), timebin = as.factor(timebin)) %>% 
  ez::ezANOVA(dv=.(scl), wid=.(vp), 
              within=.(ActPass,Intensity,timebin), 
              #between=.(pairs),
              detailed=T, type=2)


anova %>% anova_apa()
anova$ANOVA[2,] %>% partial_eta_squared_ci() #ActPass
anova$ANOVA[3,] %>% partial_eta_squared_ci() #Intensity
anova$ANOVA[4,] %>% partial_eta_squared_ci() #time
anova$ANOVA[5,] %>% partial_eta_squared_ci() #ActPass * Intensity
anova$ANOVA[6,] %>% partial_eta_squared_ci() #ActPass * Time
anova$ANOVA[7,] %>% partial_eta_squared_ci() #Intensity * time

# t-tests 
eda.mean = eda_df_long %>%
  group_by(ID, ActPass,Intensity)%>%
  summarize(value = mean(scl))
t.data.highshock= eda.mean %>%
  filter(ActPass == "Passive" & Intensity == "High")
t.data.flight = eda.mean %>%
  filter(ActPass == "Active" & Intensity == "High")
t.data.lowshock = eda.mean %>%
  filter(ActPass == "Passive" & Intensity == "Low")
t.data.avoidance = eda.mean %>%
  filter(ActPass == "Active" & Intensity == "Low")
#High shock vs low shock
t.test(t.data.highshock$value,t.data.lowshock$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)
#High shock vs. flight
t.test(t.data.highshock$value,t.data.flight$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)
#Low shock vs. avoidance
t.test(t.data.lowshock$value,t.data.avoidance$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)
#Flight vs. avoidance
t.test(t.data.avoidance$value,t.data.flight$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)


# long to wide if necessary ----------------------------------------------------

# eda_wide_out <- eda_df_long %>% 
#   group_by(ID,condition,timebin) %>%  summarise(scl = mean(scl), .groups = "drop") %>% 
#   mutate(condition = factor(condition, levels=c(1,2,3,4),labels =c("HI_threat","HI_flight","LI_threat","LI_flight")),
#          condition_names = paste0(condition,"_",timebin)) %>% select(-condition, -timebin) %>%
#   pivot_wider(., names_from = condition_names, values_from = scl)
# write.csv2(eda_wide_out,"EDA_Daten.csv",row.names=F)
# 
# #scl.data <- eda_df %>% select(ID,trial,Baseline,scl_1:scl_10)
