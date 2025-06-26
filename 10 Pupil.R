# Pupillen Analyse

#packages 
#library(signal)
library(zoo)
library(tidyverse)
library(stringr)
library(afex)

options(dplyr.summarise.inform = FALSE)

# paths
path <- paste0(getwd(),"/")
path.pupil <- paste0(path, "Data/Tobii/pupil/")
path.trigger <- paste0(path, "Data/Tobii/Trigger/")

# options & plots
cr_plots = T


sample_rate = 60  #samplerate after export
trial_length = 12  #length of trial in s

downsampling = F
sample_rate_new = 10 #for downsampling
check_ds_plots = F

standardizing = F

lowpass = T
lowPassFreq = 2 #low pass filter (Hz) entweder keinen oder 2Hz wie bei Matthias
highpass = F
highPassFreq = 0.01 #high pass filter (Hz)
check_smoothed_plots = F #not yet coded

responseWindow = c(2.5,3.5)
baselineWindow = c(-1,0)


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
}


filemat = list.files(paste0(path, "Data/Tobii/pupil/"), pattern=".*_Pupildata.txt")
filemat <- filemat[!grepl("vp05|vp22|vp30|vp51", filemat)]  #filter general exclusions
triggermat = list.files(paste0(path, "Data/Tobii/Trigger/"), pattern=".*_Trigger.dat")


pupils_list = list()
pupils_df_list = list()
pupil_df = tibble()
ga_unified = tibble()
code_times = c()


for (subject_inmat in filemat){ 
  
  #subject_inmat = filemat[[1]]
  start.time <- Sys.time()  
  
  pupil = read.csv(paste0(path.pupil,subject_inmat),dec=".", sep="\t",
                   stringsAsFactors =F, header = T) %>% #read export
    mutate(index = suppressWarnings(as.numeric(index)), pos_time = suppressWarnings(as.numeric(pos_time)),
           x_diameter = suppressWarnings(as.numeric(x_diameter)), y_diameter = suppressWarnings(as.numeric(y_diameter)),
           marker = suppressWarnings(as.numeric(marker))) %>%
    filter(!is.na(index)) %>%
    mutate(index = 1:nrow(.)) %>%
    mutate(time = (pos_time-min(pos_time))/1000, sample = 1:n()) %>% 
    rowwise() %>%  mutate(diameter = mean(x_diameter,y_diameter),
                          diameter = ifelse(diameter==-10000,NA,diameter),
                          trigger = ifelse(marker==-1,0,marker)) %>% 
    select(sample,time,diameter,trigger) 
  
  filename =  filemat[filemat == subject_inmat] %>% str_remove(.,pattern="_Pupildata.txt")# %>% str_remove(.,pattern="tfo")
  
  conditions = read.table(paste0(path.trigger,filename,"_Trigger.dat"), sep="\t",stringsAsFactors =F) %>% .$V1
  # hier dann mit trigger/SiCP_
  
  if (filename == "tfo06"){
    conditions = conditions[17:196]
  
  }
  
  
  for (i in 1:length(conditions)){pupil$trigger[pupil$trigger == i] = conditions[[i]]}
  
  # {if (sum(conditions == 1) != 40) { warning(paste0("Warning: Too few/many Threat triggers"))}
  #   if (sum(conditions == 2) != 40) { warning(paste0("Warning: Too few/many  Flight triggers"))}
  #   if (sum(conditions == 3) != 40) { warning(paste0("Warning: Too few/many  Safe triggers"))}
  #   if (sum(conditions == 4) != 40) { warning(paste0("Warning: Too few/many  Safe triggers"))}
  # }
  # 
  

  print("data reading successful!")
  
  # interpolation of missing values  
  
  pupil = pupil %>% rowwise() %>% mutate(interpolated = ifelse(is.na(diameter),T,F)) %>% ungroup() #create coloum for interpolated values
  pupil$diameter[[1]] = pupil$diameter[[min(which(!is.na(pupil$diameter)))]] # remove leading NA
  pupil$diameter[[length(pupil$diameter)]] = pupil$diameter[[max(which(!is.na(pupil$diameter)))]] # remove trailing NA
  pupil$diameter = na.approx(pupil$diameter) #interpolate NA
  
  print("interpolating missing values done!")
  
  
  # artefact correction
  
  if (downsampling) {
    
    #downsample (all columns)
    triggers_time = pupil$time[pupil$trigger != 0] #eda$time[eda$Trigger %>% is.na() == FALSE]
    conversion = round(sample_rate / sample_rate_new)
    pupil_downsampled = data.frame(time=sample.down(pupil$time,conversion),
                                   diameter=sample.down(pupil$diameter,conversion),
                                   interpolated = sample.down(pupil$interpolated,conversion),
                                   trigger=0) %>% 
      mutate(sample=1:n()) %>% select(sample, everything()) %>%
      mutate(interpolated = ifelse(interpolated < 0.5,F,T))
    
    #calculate closest position for trigger onsets in downsampled time
    triggers_time_old = pupil$time[pupil$trigger != 0]
    triggers_indices_new = triggers_time_old %>% closestRepresentative(pupil_downsampled$time, returnIndices = T) #for each old trigger time, find index of closest existing downsampled time
    pupil_downsampled$trigger[triggers_indices_new] = pupil$trigger[pupil$trigger != 0] #inject conditions as triggers of onsets
    
    # eda_downsampled %>% ggplot(aes(x=time, y=eda)) + geom_line() + geom_vline(xintercept = eda$time[eda$trigger!=0])
    if(check_ds_plots){
      pupil %>% ggplot(aes(x=time, y=diameter)) + geom_line() + geom_line(data=pupil_downsampled, color="red",linetype=4) + 
        geom_vline(xintercept = pupil$time[pupil$trigger!=0],alpha=0.1) 
    }  
    pupil = pupil_downsampled; rm(pupil_downsampled)
    print("downsampling done!")
    
  }  
  
  #smoothing / filtering
  if (lowpass) {
    reps = 100 #filters can produce artifacts on the edges => fill edges with first and last values
    diameter_filtered = c(rep(first(pupil$diameter), reps), pupil$diameter, rep(last(pupil$diameter), reps)) %>% 
      signal::filtfilt(signal::butter(2, lowPassFreq/(ifelse(downsampling, sample_rate_new, sample_rate)/2)), .) #low-pass filter
    
    if(highpass){
      diameter_filtered = diameter_filtered %>%  # high-pass filter
        signal::filtfilt(signal::butter(8, type="high", highPassFreq/(ifelse(downsampling, sample_rate_new, sample_rate)/2)), .)  
    }
    
    pupil_filt = pupil
    pupil_filt$diameter = diameter_filtered[(reps+1):(length(diameter_filtered)-reps)]; rm(diameter_filtered)
    
    
    if (check_smoothed_plots) {
      pupil  %>% ggplot(aes(x=time, y=diameter)) + geom_line() + geom_line(data=pupil_filt, color="red",linetype=4) + 
        geom_vline(xintercept = pupil$time[pupil$trigger!=0],alpha=0.1) 
    }
    
    pupil = pupil_filt
    print("low pass filtering done!")
  }  
  
  pupil_vp = pupil %>%  select(-diameter) %>%  filter(trigger != 0) %>% ungroup() %>% #One trial per row per VP
    mutate(trial = 1:n(), condition = trigger,
           time_start = time, time_end = time + trial_length,
           sample_start = sample, 
           sample_end =  {sample+(trial_length*ifelse(downsampling,sample_rate_new,sample_rate))} %>% round()) %>% 
    select(-trigger, -time, -sample, -interpolated) %>% select(trial, condition, everything())
  
  for (row in nrow(pupil_vp)){ #count interpolated samples per trial
    pupil_vp$interpolated = pupil_vp$time_start %>% lapply(function(start)
      pupil %>% filter(time >= start-abs(min(baselineWindow)), time < start + trial_length) %>% .$interpolated %>% mean(.) %>% round(.,3))
    
    exclude = rep(TRUE,ifelse(downsampling, sample_rate_new/2, sample_rate/2)) #exclude if more than 30 interpolated samplepoints (0.5s) in a row
    
    pupil_vp$valid = pupil_vp$time_start %>% lapply(function(start)
      pupil %>% filter(time >= start-abs(min(baselineWindow)), time < start + trial_length) %>% {!grepl(paste(exclude,collapse=""),paste(.$interpolated,collapse=""))})
  }
  
  
  pupils_df_list[[filename]] = pupil_vp
  pupils_list[[filename]] = pupil
  
  
  unified = pupil_vp %>% mutate(diameter = time_start %>% lapply(function(start) 
    pupil %>% filter(time >= start - abs(min(baselineWindow)),
                     time <= start + trial_length) %>% select(-trigger)))
  
  
  for (t in 1:nrow(unified)) {
    unified$diameter[[t]] = unified$diameter[[t]] %>% 
      mutate(trial = t, 
             #time = time - (min(time)+abs(min(baselineWindow))), 
             #samplepoint = sample-min(sample),
             condition = unified$condition[[t]], 
             #diameter = diameter - unified$baseline[[t]] #unify starting time to allow overlap
      )}
  
  if (standardizing) {  
    unified$diameter <- unified$diameter %>% bind_rows() %>% 
      mutate(diameter = (diameter - mean(diameter))/sd(diameter)) %>%
      #z-stand %>%
      split(.$trial)
    
    print("z-standardizing done!")
    
  }
  
  
  # baseline and level scoring  
  
  test.time <- Sys.time()  
  pupil_diameter_binded <- unified$diameter %>% bind_rows()
  
  for (t in 1:nrow(pupil_vp)) {
    start = pupil_vp$time_start[t]
    
    baseline_trial = pupil_diameter_binded %>% filter(time <= start + max(baselineWindow),
                                                                 time >= start + min(baselineWindow)) %>% 
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial = pupil_diameter_binded %>% filter(time <= start + max(responseWindow),
                                                                       time >= start + min(responseWindow)) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_1 = pupil_diameter_binded %>% filter(time >= start + 0, time <= start + 0.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_2 = pupil_diameter_binded %>% filter(time >= start + 0.5,time <= start + 1) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_3 = pupil_diameter_binded %>% filter(time >= start + 1, time <= start + 1.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_4 = pupil_diameter_binded %>% filter(time >= start + 1.5,time <= start + 2) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_5 = pupil_diameter_binded %>% filter(time >= start + 2, time <= start + 2.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_6 = pupil_diameter_binded %>% filter(time >= start + 2.5,time <= start + 3) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_7 = pupil_diameter_binded %>% filter(time >= start + 3, time <= start + 3.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_8 = pupil_diameter_binded %>% filter(time >= start + 3.5,time <= start + 4) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_9 = pupil_diameter_binded %>% filter(time >= start + 4, time <= start + 4.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_10 = pupil_diameter_binded %>% filter(time >= start + 4.5, time <= start + 5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_11 = pupil_diameter_binded %>% filter(time >= start + 5,time <= start + 5.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_12 = pupil_diameter_binded %>% filter(time >= start + 5.5, time <= start + 6) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_13 = pupil_diameter_binded %>% filter(time >= start + 6,time <= start + 6.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_14 = pupil_diameter_binded %>% filter(time >= start + 6.5, time <= start + 7) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_15 = pupil_diameter_binded %>% filter(time >= start + 7,time <= start + 7.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_16 = pupil_diameter_binded %>% filter(time >= start + 7.5, time <= start + 8) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_17 = pupil_diameter_binded %>% filter(time >= start + 8, time <= start + 8.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_18 = pupil_diameter_binded %>% filter(time >= start + 8.5, time <= start + 9) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_19 = pupil_diameter_binded %>% filter(time >= start + 9, time <= start + 9.5) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    diameter_level_trial_20 = pupil_diameter_binded %>% filter(time >= start + 9.5, time <= start + 10) %>%  
      summarize(mean = mean(diameter,na.rm=T)) %>% unlist()
    
    
    dilations = tibble(dilation = diameter_level_trial - baseline_trial,
                       dilation_1 = diameter_level_trial_1 - baseline_trial,
                       dilation_2 = diameter_level_trial_2 - baseline_trial,
                       dilation_3 = diameter_level_trial_3 - baseline_trial,
                       dilation_4 = diameter_level_trial_4 - baseline_trial,
                       dilation_5 = diameter_level_trial_5 - baseline_trial,
                       dilation_6 = diameter_level_trial_6 - baseline_trial,
                       dilation_7 = diameter_level_trial_7 - baseline_trial,
                       dilation_8 = diameter_level_trial_8 - baseline_trial,
                       dilation_9 = diameter_level_trial_9 - baseline_trial,
                       dilation_10 = diameter_level_trial_10 - baseline_trial,
                       dilation_11 = diameter_level_trial_11 - baseline_trial,
                       dilation_12 = diameter_level_trial_12 - baseline_trial,
                       dilation_13 = diameter_level_trial_13 - baseline_trial,
                       dilation_14 = diameter_level_trial_14 - baseline_trial,
                       dilation_15 = diameter_level_trial_15 - baseline_trial,
                       dilation_16 = diameter_level_trial_16 - baseline_trial,
                       dilation_17 = diameter_level_trial_17 - baseline_trial,
                       dilation_18 = diameter_level_trial_18 - baseline_trial,
                       dilation_19 = diameter_level_trial_19 - baseline_trial,
                       dilation_20 = diameter_level_trial_20 - baseline_trial,
                       baseline = baseline_trial,
                       mean_diameter = diameter_level_trial) %>% select(dilation:dilation_20, baseline, mean_diameter, everything()) 
    
    
    pupil_vp[t, names(dilations)] = dilations
  } 
  end.time <- Sys.time() - test.time
  
  print("0.5s bin averaging done!")
  print(end.time)
  
  pupil_vp = pupil_vp %>% group_by(condition) %>% mutate(trial_condition = row_number()) %>% ungroup()
  pupils_df_list[[filename]] = pupil_vp
  
  
  
  for (t in 1:nrow(unified)) {
    unified$diameter[[t]] = unified$diameter[[t]] %>% 
      mutate(time = time - (min(time)+abs(min(baselineWindow))), 
             samplepoint = sample-min(sample),
             diameter = diameter - pupil_vp$baseline[[t]] #unify starting time to allow overlap
      )}
  
  if (cr_plots) {
    
    
    # unified$diameter %>% bind_rows() %>% 
    #   mutate(trial = as.factor(trial)) %>% 
    #   {ggplot(., aes(x=time, y=diameter, color=trial)) + facet_wrap(vars(condition)) +
    #       geom_vline(xintercept=2, color="blue",linetype="dashed") + #borders of min/max scoring
    #       geom_vline(xintercept=0, color="black",linetype="solid") + #zero 
    #       geom_path() + scale_color_viridis_d() +
    #       ggtitle(filename) + theme_bw() + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))} %>% 
    #   #print()
    #   # ggsave(paste0("E:/Experimente/Experiment - TFO/Plots/Pupille/single CS/", filename, ".png"), plot=., 
    #   #        device="png", width=1920/150, height=1080/150, dpi=300)
    # 
    # 
    # 
    # unified$diameter %>% bind_rows() %>% 
    #   mutate(condition = as.factor(condition)) %>% 
    #   group_by(condition,samplepoint) %>% 
    #   summarise(diameter = mean(diameter), time = mean(time)) %>%
    #   {ggplot(., aes(x=time, y=diameter, color=condition, group=condition,linetype=condition)) +
    #       geom_vline(xintercept=2, color="blue",linetype="dashed") + #borders of min/max scoring
    #       geom_vline(xintercept=0, color="black",linetype="solid") + #zero 
    #       geom_path() + scale_color_viridis_d() +
    #       ggtitle(filename) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))} %>% 
    #   #print()
    #   # ggsave(paste0("E:/Experimente/Experiment - TFO/Plots/Pupille//mean CS/", filename, ".png"), plot=., 
    #   #        device="png", width=1920/150, height=1080/150, dpi=300)
    # 
    # print("Cr plotting successful!")
    # 
  }
  
  unified = unified %>% mutate(ID = filename)
  ga_unified = rbind(ga_unified,unified)
  pupil_df  = rbind(pupil_df,pupil_vp %>% mutate(ID=filename))
  
  
  # stuff for calculating computing duration
  
  end.time <- Sys.time()  
  time.taken <- end.time - start.time
  code_times = c(code_times,round(as.numeric(time.taken),1))
  mean_code_time = mean(code_times)
  
  
  if(match(subject_inmat,filemat) == 1){
    print(paste0(round(match(subject_inmat,filemat)/length(filemat)*100,1), "% ..... ",
                 match(subject_inmat,filemat), " of ", length(filemat), " files processed!"))}
  else{
    print(paste0(round(match(subject_inmat,filemat)/length(filemat)*100,1), "% ..... ",
                 match(subject_inmat,filemat), " of ", length(filemat), " files processed! ..... ", 
                 round((mean_code_time*length(filemat)-mean_code_time * match(subject_inmat,filemat))/60,1), 
                 " min remaining"))
  }
  
} # end inmat for loop


pupil_df <- pupil_df %>%
  mutate(ID = gsub("_1","", ID))%>%
  mutate(ID = gsub("_2","", ID))




# Read or Save ga_unified and pupil_df for further processing

# saveRDS(ga_unified,"ET_ga_unified.RData")
# saveRDS(pupil_df,"ET_pupil_df.RData")

#######################################################
## FOR DATA ANALYSIS READ PREPROCESSED FILES HERE #####
#######################################################

ga_unified <- readRDS("ET_ga_unified.RData")
pupil_df <- readRDS("ET_pupil_df.Rdata")

# Signal quality check

mean(as.numeric(pupil_df$interpolated)) #30.4% interpolated sample points
sum(pupil_df$valid == T) #8544 valid trials
sum(pupil_df$valid == F) #4536 invalid trials
sum(pupil_df$valid == T) / (sum(pupil_df$valid == T) + sum(pupil_df$valid == F) ) # 65.3% valid trials

pupil_df %>% group_by(condition) %>%
  summarise(valid_sp = mean(as.numeric(interpolated)), 
            valid_trials = sum(valid == T))  # no condition differences!

# Plot Unified Data

# ga_unified %>% filter(valid == TRUE) %>%
#   #filter(ID %in% responder) %>%
#   .$diameter %>% bind_rows() %>%
#   mutate(condition = as.factor(condition)) %>% 
#  # filter(trial < 60) %>%
#   filter(condition %in% c(1,2,3,4)) %>% 
#   group_by(condition,samplepoint) %>% 
#   summarise(diameter = mean(diameter), time = mean(time)) %>%
#   {ggplot(., aes(x=time, y=diameter, color=condition, group=condition)) +
#       #geom_vline(xintercept=c(responseWindow), color="blue",linetype="dashed") + #borders of min/max scoring
#       #geom_vline(xintercept=0, color="black",linetype="solid") + #zero 
#       #geom_vline(xintercept=2, color="lightblue",linetype="dashed") + #US time
#       #geom_vline(xintercept=10, color="lightblue",linetype="dashed") + #US time
#       geom_path() + 
#       scale_x_continuous("Time [s]",limits=c(-1, 10)) +
#       scale_color_manual(values=c("red","orange","darkgreen", "lightgreen"), labels = c("HI Shock","HI Flight","LI Shock", "LI Flight")) +
#       scale_y_continuous("Pupil Diameter", limits = c(-.50, .25)) +
#       theme_classic() +
#       theme(legend.position=c(0.15,0.88),
#             legend.title = element_blank(),
#             #panel.grid.major.y = element_line(size=0.5, color="#DDDDDD")
#       ) }
#ggsave("../Plots/Pupille/cs_acq.png",type="cairo-png", width=1920/400, height=1080/300, dpi=300)


# Inference statistics--------------------------------------------------------------------

# read in prot files 
All_Prot_Data <- read.csv2(paste0(path.prot,"All_Prots.csv"))

# filter out problematic trials 

pupil <- pupil_df %>% rename(vp = ID) %>% filter(vp != "vp20") %>% group_by(vp)%>% mutate(trial = 1:160)
pupil_vp20 <- pupil_df %>% rename(vp = ID) %>% filter(vp == "vp20") %>% group_by(vp)%>% mutate(trial = 1:140)
pupil <- rbind(pupil,pupil_vp20)
pupil <- full_join(All_Prot_Data, pupil)

pupil <- pupil %>%
  filter(problem == 0) %>%
  mutate(condition = as.factor(condition))%>%
  mutate(ActPass = case_when(condition ==1 | condition == 3 ~ "Passive",
                             condition == 2 | condition == 4 ~ "Active"),
         Intensity = case_when(condition == 1 | condition == 2 ~ "High",
                               condition== 3 | condition == 4 ~ "Low")) 

pupil_long <- pupil %>%
  select(vp,dilation_1:dilation_20,cue,trial_condition) %>%
  pivot_longer(dilation_1:dilation_20, names_to = "bin", values_to ="value") %>%
  mutate(bin = gsub("dilation_", "", bin))%>%
  mutate(bin = as.numeric(bin)/2)%>% 
  select(-trial_condition)%>%
  group_by(vp, cue, bin) %>%
  summarize(value = mean(value, na.rm =T))%>%
  mutate(cue = as.factor(cue))%>%
  mutate(ActPass = case_when(cue ==1 | cue == 3 ~ "Passive",
                             cue == 2 | cue == 4 ~ "Active"),
         Intensity = case_when(cue == 1 | cue == 2 ~ "High",
                               cue== 3 | cue == 4 ~ "Low")) 
names(pupil_long) <- c("code","cue", "bin","value", "ActPass", "Intensity")

write.csv2(pupil_long, paste0(path.pupil, "pupil_long.csv"))


#ANOVA
anova_pupil <- pupil_long %>%
  filter(bin > 2) %>%
  mutate(code = as.factor(code), ActPass = as.factor(ActPass), Intensity = as.factor(Intensity), bin = as.factor(bin)) %>%
  ez::ezANOVA(dv=.(value), wid=.(code), 
              within=.(ActPass, Intensity, bin), 
              #between=.(pairs),
              detailed=T, type=2)
  
anova_pupil %>% apa::anova_apa()

anova_pupil$ANOVA[2,] %>% partial_eta_squared_ci() #ActPass
anova_pupil$ANOVA[3,] %>% partial_eta_squared_ci() #Intensity
anova_pupil$ANOVA[4,] %>% partial_eta_squared_ci() #bin
anova_pupil$ANOVA[5,] %>% partial_eta_squared_ci() #ActPass * Intensity
anova_pupil$ANOVA[6,] %>% partial_eta_squared_ci() #ActPass * bin
anova_pupil$ANOVA[7,] %>% partial_eta_squared_ci() #Intensity * bin
anova_pupil$ANOVA[8,] %>% partial_eta_squared_ci() #ActPass * Intensity * bin

# t-tests 
pupil.mean = pupil_long %>%
  group_by(code, ActPass,Intensity)%>%
  summarize(value = mean(value))
t.data.highshock= pupil.mean %>%
  filter(ActPass == "Passive" & Intensity == "High")
t.data.flight = pupil.mean %>%
  filter(ActPass == "Active" & Intensity == "High")
t.data.lowshock = pupil.mean %>%
  filter(ActPass == "Passive" & Intensity == "Low")
t.data.avoidance = pupil.mean %>%
  filter(ActPass == "Active" & Intensity == "Low")
#High shock vs low shock
t.test(t.data.highshock$value,t.data.lowshock$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)
#High shock vs. flight
t.test(t.data.highshock$value,t.data.flight$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)
#Low shock vs. avoidance
t.test(t.data.lowshock$value,t.data.avoidance$value,alternative = c("less"), paired=T) %>% apa::t_apa(es_ci=T)
#Flight vs. avoidance
t.test(t.data.avoidance$value,t.data.flight$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)

# long to wide if necessary ----------------------------------------------------

pupil_wide_out <- pupil_df_long %>%
  #filter(ID %in% responder) %>%
  select(ID, condition, timebin, diameter) %>%
  group_by(ID,condition,timebin) %>%  summarise(diameter = mean(diameter, na.rm = T), .groups = "drop") %>%
  mutate(condition = factor(condition, levels=c(1,2,3,4),labels =c("HI_threat","HI_flight","LI_threat", "LI_flight")),
         condition_names = paste0(condition,"_",timebin)) %>% select(-condition, -timebin) %>%
  pivot_wider(., names_from = condition_names, values_from = diameter)

write.csv2(pupil_wide_out,"Pupillen_Daten.csv",row.names=F)




