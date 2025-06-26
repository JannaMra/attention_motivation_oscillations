############################################################################
# ShockScanning_Flight Oscillations project
# 
# Get Fixations from Gazedata-Files
# Design:
# cue 1: high intensity shock (shock)
# cue 2: high intensity flight (flight)
# cue 3: low intensity shock (safety)
# cue 4: low intensity flight (active control)
# Stimulation:
# 2s fixation cross -> 8s stimulus -> 6-8s fixation cross with shock or not

library(tidyverse)
library(readr)

path <- "C:/Users/jat41gk/Documents/Projekte/Action-Motivation-Oscillations/Analyse/" #paste0(getwd(),"/")
path.log <- "Data/log/" %>% paste0(path, .) 
path.seq <- "Data/sequences/" %>% paste0(path, .) 
path.eye <- "Data/Tobii/eye/" %>% paste0(path, .) 
path.trigger <- "Data/Tobii/Trigger/" %>% paste0(path, .) 

files.eye = list.files(path.eye, pattern=".*.txt")
files.trigger = list.files(path.trigger, pattern=".*.dat")
  

{ # Global Parameters -------------------------------------------------------
  exclusions.eye.num = c() #a priori exclusions, e.g. calibration not successful
  numToEye = function(nums) return(nums %>% formatC(width=2, format="d", flag="0") %>% paste0("vp", .)) #add leading zeros and prefix "vp"
  exclusions.eye = exclusions.eye.num %>% numToEye()
  
  #raw data parsing
  saccadeAmpMin = 1 #minimum saccade amplitude in pixels (should correspond to .5 - 1 degree of visual angle); smaller saccades will lead to a combination of adjacent fixations
  parsingTrialPlots = F
  
  pathToCode = function(path, path.sep="/", file.ext="\\.") {
    first = path %>% gregexpr(path.sep, .) %>% sapply(max) %>% {. + 1}
    last = path %>% gregexpr(file.ext, .) %>% sapply(max) %>% {. - 1}
    return(path %>% substring(first, last))
  }
  
  trials.n = 160 #number of trials that shall be analyzed (if more trials, last ones will be taken)
  
  preStim = 2000
  trialEnd = 10000
  
  screen.height = 1024 #height of screen in pix
  screen.width  = 1280 # width of screen in pix
}


# Read & parse data -------------------------------------------------------
fixations = data.frame()
messages = data.frame()
setwd(path.eye)
for (filename in files.eye) {
  file = read_table(filename) #read.delim(filename)# read_table
  subjectName = pathToCode(filename) %>% gsub("_Gazedata","", .)
  cat(subjectName, "\n")
  messages.s = file %>% filter(marker != -1) %>% transmute(subject = subjectName, trial = marker, pos_time = pos_time)
  trialStarts = messages.s %>% .$pos_time
  
  fixations.s = data.frame()
  for (trial in seq_along(trialStarts)) {
    trialStart = trialStarts[trial]
    file.trial = file %>% filter(pos_time >= trialStart-500 , pos_time <= trialStart + trialEnd) %>% 
      mutate(pos_x = ifelse(pos_x == -10000, NA, pos_x), pos_y = ifelse(pos_y == -10000, NA, pos_y))
    blinks = file.trial %>% with(pos_x %>% is.na() | pos_y %>% is.na()) %>% as.integer()
    fixations.trial = file.trial %>% with(eventdetection_lowspeed(pos_time, pos_x, pos_y, blinks, saccadeAmpMin, plotdata=parsingTrialPlots)) %>% 
      .$fixations %>% data.frame() %>% mutate(TRIAL_LABEL = trial) %>% select(TRIAL_LABEL, everything()) #%>% rename(trial = TRIAL_LABEL, start = X1, end = X2, x = X3, y = X4)
    
    if (fixations.trial %>% ncol() == 5) {
      names(fixations.trial) = c("trial","start","end","x","y")
      fixations.s = fixations.s %>% bind_rows(fixations.trial)
    } else {
      #fixations.s = fixations.s %>% bind_rows(data.frame(trial=trial, start=NA, end=NA, x=NA, y=NA))
    }
    if (parsingTrialPlots) readline(paste0(subjectName, ", trial ", trial, " (press enter to continue)"))
  }
  
  fixations.s = fixations.s %>% mutate(subject = subjectName) %>% select(subject, everything())
  fixations = fixations %>% bind_rows(fixations.s)
  messages = messages %>% bind_rows(messages.s)

}

fixations <- fixations %>%
  mutate(block = ifelse(grepl("_2", fixations$subject, fixed = TRUE), 2, 1))%>%
  mutate(trial = ifelse(grepl("_2", fixations$subject, fixed = TRUE), trial+80, trial))%>%
  mutate(x = x + screen.width/2, #center of screen = (0, 0) => switch to bottom left = (0, 0)
         y = y + screen.height/2)
fixations$subject <- gsub("\\_[1-2]","",fixations$subject)
messages <- messages %>%
  mutate(block = ifelse(grepl("_2", messages$subject, fixed = TRUE), 2, 1))%>%
  mutate(trial = ifelse(grepl("_2", messages$subject, fixed = TRUE), trial+80, trial))%>%
  rename(time = pos_time)
messages$subject <- gsub("\\_[1-2]","",messages$subject)
  
# messages = messages %>% 
#   rename(time = pos_time) %>% #for validating baselines
#   mutate(subject = subject %>% sub("vp", "", .) %>% as.integer(),
#          #block = block %>% as.numeric(),
#          #trial = trial %>% sub("Trial: ","", .) %>% as.numeric(), #for EyeLink
#          trial = trial)

# fixations = fixations %>% 
#   mutate(subject = subject %>% sub("vp", "", .) %>% as.integer(),
#          x = x + screen.width/2, #center of screen = (0, 0) => switch to bottom left = (0, 0)
#          y = y + screen.height/2,
#          #y = screen.height - y, #revert y-axis
#          #trial = trial %>% sub("Trial: ","", .) %>% as.numeric(),
#          trial = trial)


# Write tab delimiter text file
write.table(fixations, "Data/Tobii/Fixations.txt" %>% paste0(path, .)) 
write.table(messages, "Data/Tobii/Messages.txt" %>% paste0(path, .))



