############################################################################
# Analyze ShockScanning_Flight project
# 
# Calculate fixation numbers and durations as a function of experimental conditions
# Design:
# cue 1: high intensity shock (shock)
# cue 2: high intensity flight (flight)
# cue 3: low intensity shock (safety)
# cue 4: low intensity flight (active control)
# Stimulation:
# 2s fixation cross -> 8s stimulus -> 6-8s fixation cross with shock or not

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# General: loading packages, definition of variables, subjects etc.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
require("lme4")
require("lmerTest")
require("performance")
require("BayesFactor")

#rm(list=ls()) 

# paths
path <- paste0(getwd(),"/")

# Onsets laden
msg <- read.table(paste(path,"Data/Tobii/Messages.txt",sep=""))
names(msg) <- c("vp","trial","time")
fixa <- read.table(paste(path,"Data/Tobii/Fixations.txt", sep=""))

# Exclusions:
exclusions <- c("vp05","vp22","vp30","vp51")
eye.invalid.bl <- c("vp04", "vp08", "vp09", "vp20", "vp31", "vp46")

# Determine which subjects should be analyzed
vpn = fixa$subject %>% unique() %>% sort() %>% as.character() #all subjects in fixations
vpn.n = length(vpn)
vpn = vpn[!(vpn %in% c("vp22"))] #vp22 not enough values per condition
vpn <-  vpn[!(vpn %in% exclusions)]
vpn <-  vpn[!(vpn %in% eye.invalid.bl)]
vpn.n = length(vpn)

All_Data <- data.frame()

# combine data in one dataset---------------------------------------------------

All_Prot_Data <- read.csv2(paste0(path, "Data/prot/All_Prots.csv"))%>%
  rename(subject = vp)

# Eyetracking data
All_Eye_Data <- read.csv2(paste0(path, "Data/Tobii/All_Eye_Data.csv"))
All_Eye_Data <- All_Eye_Data %>%
  mutate(cb = rowMeans(All_Eye_Data%>% select(X5_CB:X8_CB), na.rm=TRUE),
         fixnum = rowMeans(All_Eye_Data%>% select(X5_FixN:X8_FixN), na.rm=TRUE),
         fixdur = rowMeans(All_Eye_Data%>% select(X5_FixDur:X8_FixDur), na.rm=TRUE))%>%
  rename(subject = vp)%>%
  select(subject, trial,cue, rt, problem, cb, fixnum, fixdur)

# add HR data
hr <- read_rds("heart_df.rds") 
hr$time <- as.numeric(unlist(hr$time))
hr <- hr %>%
  filter(time > 6 & time < 11)%>%
  group_by(subject, trial)%>%
  summarize(hr = mean(HRchange))%>%
  ungroup()%>%
  mutate(subject = paste0("vp", ifelse(subject < 10, "0",""),subject))

# add SCL data
eda <-  readRDS("EDA_unified.RData")
eda <- eda %>%
  rename(subject = ID)%>%
  mutate(scl = rowMeans(eda%>% select(scl_7:scl_10), na.rm=TRUE))%>%
  select(subject, trial, scl)

# add pupil data
pupil <-  readRDS("ET_pupil_df.RData")%>%
  rename(subject = ID)%>%
  group_by(subject)%>%
  mutate(trial = 1:length(trial))%>%
  ungroup()

pupil <- pupil %>%
  mutate(dilation = rowMeans(pupil %>% select(dilation_13:dilation_20), na.rm=TRUE))%>%
  select(subject, trial, dilation)

# add EEG data
eeg <- read.csv2("./Data/EEG/single_trial_alpha_all_conditions.csv")%>%
  rename(subject = ID)%>%
  select(subject, trial, alpha_bl)%>%
  rename(alpha = alpha_bl)

# add high alpha data
eeg_high <- read.csv2("./Data/EEG/single_trial_high_alpha_all_conditions.csv")%>%
  rename(subject = ID)%>%
  select(subject, trial, alpha_bl)%>%
  rename(alpha_high = alpha_bl)

# add low alpha data
eeg_low <- read.csv2("./Data/EEG/single_trial_low_alpha_all_conditions.csv")%>%
  rename(subject = ID)%>%
  select(subject, trial, alpha_bl)%>%
  rename(alpha_low = alpha_bl)


All_Data <- full_join(All_Prot_Data, All_Eye_Data) %>% 
  full_join(hr, c("subject", "trial")) %>% 
  full_join(eda, by=c("subject", "trial")) %>% 
  full_join(pupil, by=c("subject", "trial")) %>%
  full_join(eeg, by = c("subject", "trial")) %>%
  full_join(eeg_high, by = c("subject", "trial")) %>%
  full_join(eeg_low, by = c("subject", "trial")) %>%
  filter(!subject %in% exclusions)%>%
  filter(subject != "vp22") #vp22 not enough values per condition

# LMM RTs ----------------------------------------------------------------------
All_Data <- All_Data %>%
  filter(problem == 0)%>%
  mutate(cb = datawizard::standardize(cb), fixnum = datawizard::standardize(fixnum), fixdur = datawizard::standardize(fixdur),hr = datawizard::standardize(hr), dilation = datawizard::standardize(dilation), alpha = datawizard::standardize(alpha), alpha_high = datawizard::standardize(alpha_high), alpha_low = datawizard::standardize(alpha_low),scl = datawizard::standardize(scl))%>% #z-transformation
  #mutate(cb = scale(cb), fixdur = scale(fixdur), fixnum = scale(fixnum), hr = scale(hr), dilation = scale(dilation), alpha = scale(alpha_bl))%>% #z-transformation
  filter(cue == 2 | cue == 4)%>%
  filter(rt > 100 & rt < 1000)%>%
  filter(alpha >= -4)

model_full_simple <- lmer(rt ~ cb + hr + scl + dilation + alpha + cue + trial + (1|subject), All_Data) %>% summary() %>% print() #all

#performance::r2(model_full)
lmer(rt ~ cb + hr + scl + dilation + alpha  + (1|subject), All_Data %>% filter(cue == 2))%>% summary() %>% print() #flight painful
lmer(rt ~ cb + hr + scl + dilation + alpha  + (1|subject), All_Data %>% filter(cue == 4))%>% summary() %>% print() #flight non-painful

lmer(rt ~ cb + (1|subject), All_Data)%>% summary() %>% print()
lmer(rt ~ hr + (1|subject), All_Data)%>% summary() %>% print()
lmer(rt ~ fixnum + (1|subject), All_Data)%>% summary() %>% print()
lmer(rt ~ fixdur + (1|subject), All_Data)%>% summary() %>% print()
lmer(rt ~ dilation + (1|subject), All_Data)%>% summary() %>% print()
lmer(rt ~ scl + (1|subject), All_Data)%>% summary() %>% print()
lmer(rt ~ alpha + (1|subject), All_Data)%>% summary() %>% print()
lmer(rt ~ alpha_high + (1|subject), All_Data)%>% summary() %>% print()
lmer(rt ~ alpha_low + (1|subject), All_Data)%>% summary() %>% print()

# Exploratory ##################################################################

#Flight Painful
lmer(rt ~ cb + (1|subject), All_Data %>% filter(cue == 2)) %>% summary() %>% print()
lmer(rt ~ hr + (1|subject), All_Data %>% filter(cue == 2))%>% summary() %>% print()
lmer(rt ~ fixnum + (1|subject), All_Data %>% filter(cue == 2))%>% summary() %>% print()
lmer(rt ~ fixdur + (1|subject), All_Data %>% filter(cue == 2))%>% summary() %>% print()
lmer(rt ~ dilation + (1|subject), All_Data %>% filter(cue == 2))%>% summary() %>% print()
lmer(rt ~ scl + (1|subject), All_Data %>% filter(cue == 2))%>% summary() %>% print()
lmer(rt ~ alpha + (1|subject), All_Data %>% filter(cue == 2))%>% summary() %>% print()

#Flight Non-Painful
lmer(rt ~ cb + (1|subject), All_Data %>% filter(cue == 4))%>% summary() %>% print()
lmer(rt ~ hr + (1|subject), All_Data %>% filter(cue == 4))%>% summary() %>% print()
lmer(rt ~ fixnum + (1|subject), All_Data %>% filter(cue == 4))%>% summary() %>% print()
lmer(rt ~ fixdur + (1|subject), All_Data %>% filter(cue == 4))%>% summary() %>% print()
lmer(rt ~ dilation + (1|subject), All_Data %>% filter(cue == 4))%>% summary() %>% print()
lmer(rt ~ scl + (1|subject), All_Data %>% filter(cue == 4))%>% summary() %>% print()
lmer(rt ~ alpha + (1|subject), All_Data %>% filter(cue == 4))%>% summary() %>% print()

#Cue
lmer(rt ~ cb*cue + (1|subject), All_Data) %>% summary() %>% print()
lmer(rt ~ hr*cue + (1|subject), All_Data) %>% summary() %>% print()
lmer(rt ~ fixnum*cue + (1|subject), All_Data)%>% summary() %>% print()
lmer(rt ~ fixdur*cue + (1|subject), All_Data)%>% summary() %>% print()
lmer(rt ~ dilation*cue + (1|subject), All_Data)%>% summary() %>% print()
lmer(rt ~ scl*cue + (1|subject), All_Data)%>% summary() %>% print()
lmer(rt ~ alpha*cue + (1|subject), All_Data)%>% summary() %>% print()
