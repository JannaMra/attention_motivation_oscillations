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

## Bayesian Approach ------------------------------------------------------------

library(tidyverse)
library("scales")
library(scales)
library(patchwork)
requirePackage("BayesFactor")

# EEG Data ---------------------------------------------------------------------

alpha_long <- read.csv2("alpha_long.csv")%>%
  group_by(ID, time, ActPass, Intensity)%>%
  summarize(
    value = mean(value)
  )%>%
  mutate(
    ID = as.factor(ID),
    cue = case_when(ActPass == "Passive" & Intensity == "High" ~ 1,
                    ActPass == "Active" & Intensity == "High" ~ 2,
                    ActPass == "Passive" & Intensity == "Low" ~ 3,
                    ActPass == "Active" & Intensity == "Low" ~ 4),
    predictor_exclusive_controllability = case_when(cue == 1 ~ -1,
                                                    cue == 2 ~ 1,
                                                    cue == 3 ~ -1,
                                                    cue == 4 ~ 1),
    predictor_exclusive_intensity = case_when(cue == 1 ~ 1,
                                              cue == 2 ~ 1,
                                              cue == 3 ~ -1,
                                              cue == 4 ~ -1),
    predictor_additive = case_when(cue == 1 ~ 0,
                                   cue == 2 ~ 2,
                                   cue == 3 ~ -2,
                                   cue == 4 ~ 0),
    predictor_weighted_add = case_when(cue == 1 ~ -1,
                                       cue == 2 ~ 3,
                                       cue == 3 ~ -3,
                                       cue == 4 ~ 1),
    predictor_interactive = case_when(cue == 1 ~ -3, #-1 #-2
                                      cue == 2 ~ 5,  #4  #4
                                      cue == 3 ~ -4, #-3 #-3
                                      cue == 4 ~ 2)) #0  #2



# Prepare Dataframe for Bayes Factors
alpha_models <- data.frame(matrix(ncol = 21, nrow = 4))
names <- c("model", paste0(1:20))
colnames(alpha_models) <- names
alpha_models$model <- c("controllability", "intensity", "interactive", "controllability_plus_interaction") #"additive", "weighted", 

for (timebin in 1:20) {
  
  # Model 1: Controllability
  set.seed(30)
  BFfull_exclusive_controllability <- generalTestBF(value ~ predictor_exclusive_controllability + ID, whichRandom="ID", data=alpha_long%>% filter(time == timebin), 
                                                    neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability <- BFfull_exclusive_controllability[1]/BFfull_exclusive_controllability[2]
  extractBF(BF_controllability)$bf
  alpha_models[1,timebin+1] <- extractBF(BF_controllability)$bf
  
  # Model 2: Intensity
  set.seed(30)
  BFfull_exclusive_intensity <- generalTestBF(value ~ predictor_exclusive_intensity + ID, whichRandom="ID", data=alpha_long%>% filter(time == timebin), 
                                              neverExclude="ID", iterations=100000, whichModels = "all")
  BF_intensity <- BFfull_exclusive_intensity[1]/BFfull_exclusive_intensity[2]
  alpha_models[2,timebin+1] <- extractBF(BF_intensity)$bf
  
  # # Model 3: Controllability + Intensity
  # set.seed(30)
  # BFfull_additive <- generalTestBF(value ~ predictor_additive + ID, whichRandom="ID", data= alpha_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_additive <- BFfull_additive[1]/BFfull_additive[2]
  # 
  # #BF_additive_vs_controllability <- BFfull_additive[1]/BFfull_exclusive_controllability[1]
  # alpha_models[3,timebin+1] <- extractBF(BF_additive)$bf
  
  # # Model 4: Controllability + Intensity Weighted
  # set.seed(30)
  # BFfull_weighted <- generalTestBF(value ~ predictor_weighted_add + ID, whichRandom="ID", data=alpha_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_weighted <- BFfull_weighted[1]/BFfull_weighted[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # alpha_models[4,timebin+1] <- extractBF(BF_weighted)$bf
  
  # # Model 5: Controllability * Intensity Interaction
  # set.seed(30)
  # BFfull_interactive <- generalTestBF(value ~ predictor_interactive + ID, whichRandom="ID", data=alpha_long%>% filter(time == timebin), 
  #                                     neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_interactive <- BFfull_interactive[1]/BFfull_interactive[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # alpha_models[3,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 6: Interaction Effect 
  set.seed(30)
  BFfull_interactive <- generalTestBF(value ~ predictor_exclusive_controllability * predictor_exclusive_intensity + ID, whichRandom="ID", data=alpha_long%>% filter(time == timebin), 
                                      neverExclude="ID", iterations=100000, whichModels = "all")
  BF_interactive <- BFfull_interactive[7]/BFfull_interactive[8]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  alpha_models[3,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 7: Controllability plus Interaction Effect 
  set.seed(30)
  BFfull_controllability_interaction <- generalTestBF(value ~ predictor_exclusive_controllability + predictor_exclusive_controllability:predictor_exclusive_intensity + ID, whichRandom="ID", data=alpha_long%>% filter(time == timebin), 
                                          neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability_interaction <- BFfull_controllability_interaction[3]/BFfull_controllability_interaction[4]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  alpha_models[4,timebin+1] <- extractBF(BF_controllability_interaction)$bf
}

#alpha_models <- alpha_models %>%
#  rbind(c("interactive_vs_controllability", as.numeric(alpha_models[5, 2:21] / alpha_models[1, 2:21])))

alpha_models_long <- alpha_models %>%
  #filter(model != "interactive_vs_controllability")%>%
  pivot_longer(
    !model,
    values_to = "BF",
    names_to = "timebin"
  )%>%
  mutate(
    time = as.numeric(timebin)/2-0.25,
    BF=as.numeric(BF)
  )

alpha_models_long$model <- factor(alpha_models_long$model, levels = c("controllability", "intensity", "interactive", "controllability_plus_interaction")) #"additive", "weighted", 

bayes_alpha_plot <- ggplot(alpha_models_long, aes(x = time, y= BF,  group = model, linetype = model, color = model)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE", color = NA)+
  geom_hline(yintercept=1, color = "black")+ 
  geom_line(linewidth = 0.8)+
  scale_linetype_manual(name="Model",  values=c("solid", "solid", "dotdash", "dotdash"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model"))+
  scale_color_manual(name="Model", values = c("#239BE7","#3434AD","#38d4a2", "#1c5c5a"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model")) +
  scale_y_continuous(name = "Bayes Factor", trans='log2')+  #breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000)
  scale_x_continuous(name = "Time",limits = c(0.2,10), breaks = seq(0,10,2))+
  theme_classic()

# HR Data ----------------------------------------------------------------------

hr_long <- read_rds("heart_df.rds") %>%
  rename(ID = vp)%>%
  group_by(ID, time, ActPass, Intensity)%>%
  summarize(
    value = mean(HRchange, na.rm = T)
  )%>%
  mutate(
    ID = as.factor(ID),
    cue = case_when(ActPass == "Passive" & Intensity == "High" ~ 1,
                    ActPass == "Active" & Intensity == "High" ~ 2,
                    ActPass == "Passive" & Intensity == "Low" ~ 3,
                    ActPass == "Active" & Intensity == "Low" ~ 4),
    predictor_exclusive_controllability = case_when(cue == 1 ~ -1,
                                                    cue == 2 ~ 1,
                                                    cue == 3 ~ -1,
                                                    cue == 4 ~ 1),
    predictor_exclusive_intensity = case_when(cue == 1 ~ 1,
                                              cue == 2 ~ 1,
                                              cue == 3 ~ -1,
                                              cue == 4 ~ -1),
    predictor_additive = case_when(cue == 1 ~ 0,
                                   cue == 2 ~ 2,
                                   cue == 3 ~ -2,
                                   cue == 4 ~ 0),
    predictor_weighted_add = case_when(cue == 1 ~ -1,
                                       cue == 2 ~ 3,
                                       cue == 3 ~ -3,
                                       cue == 4 ~ 1),
    predictor_interactive = case_when(cue == 1 ~ -3, #-1 #-2
                                      cue == 2 ~ 5,  #4  #4
                                      cue == 3 ~ -4, #-3 #-3
                                      cue == 4 ~ 2)) #0  #2

# Prepare Dataframe for Bayes Factors
hr_models <- data.frame(matrix(ncol = 16, nrow = 4))
names <- c("model", paste0(1:15))
colnames(hr_models) <- names
hr_models$model <- c("controllability", "intensity", "interactive", "controllability_plus_interaction") #"additive", "weighted",

for (timebin in 1:15) {
  
  # Model 1: Controllability
  set.seed(30)
  BFfull_exclusive_controllability <- generalTestBF(value ~ predictor_exclusive_controllability + ID, whichRandom="ID", data=hr_long%>% filter(time == timebin), 
                                                    neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability <- BFfull_exclusive_controllability[1]/BFfull_exclusive_controllability[2]
  extractBF(BF_controllability)$bf
  hr_models[1,timebin+1] <- extractBF(BF_controllability)$bf
  
  # Model 2: Intensity
  set.seed(30)
  BFfull_exclusive_intensity <- generalTestBF(value ~ predictor_exclusive_intensity + ID, whichRandom="ID", data=hr_long%>% filter(time == timebin), 
                                              neverExclude="ID", iterations=100000, whichModels = "all")
  BF_intensity <- BFfull_exclusive_intensity[1]/BFfull_exclusive_intensity[2]
  hr_models[2,timebin+1] <- extractBF(BF_intensity)$bf
  
  # # Model 3: Controllability + Intensity
  # set.seed(30)
  # BFfull_additive <- generalTestBF(value ~ predictor_additive + ID, whichRandom="ID", data= hr_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_additive <- BFfull_additive[1]/BFfull_additive[2]
  # 
  # #BF_additive_vs_controllability <- BFfull_additive[1]/BFfull_exclusive_controllability[1]
  # hr_models[3,timebin+1] <- extractBF(BF_additive)$bf
  # 
  # # Model 3: Controllability + Intensity Weighted
  # set.seed(30)
  # BFfull_weighted <- generalTestBF(value ~ predictor_weighted_add + ID, whichRandom="ID", data=hr_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_weighted <- BFfull_weighted[1]/BFfull_weighted[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # hr_models[4,timebin+1] <- extractBF(BF_weighted)$bf
  # 
  # # Model 6: Controllability * Intensity Interaction
  # set.seed(30)
  # BFfull_interactive <- generalTestBF(value ~ predictor_interactive + ID, whichRandom="ID", data=hr_long%>% filter(time == timebin), 
  #                                     neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_interactive <- BFfull_interactive[1]/BFfull_interactive[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # hr_models[3,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 6: Interaction Effect 
  set.seed(30)
  BFfull_interactive <- generalTestBF(value ~ predictor_exclusive_controllability * predictor_exclusive_intensity + ID, whichRandom="ID", data=hr_long%>% filter(time == timebin), 
                                      neverExclude="ID", iterations=100000, whichModels = "all")
  BF_interactive <- BFfull_interactive[7]/BFfull_interactive[8]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  hr_models[3,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 7: Controllability plus Interaction Effect 
  set.seed(30)
  BFfull_controllability_interaction <- generalTestBF(value ~ predictor_exclusive_controllability + predictor_exclusive_controllability:predictor_exclusive_intensity + ID, whichRandom="ID", data=hr_long%>% filter(time == timebin), 
                                                      neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability_interaction <- BFfull_controllability_interaction[3]/BFfull_controllability_interaction[4]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  hr_models[4,timebin+1] <- extractBF(BF_controllability_interaction)$bf
}


#hr_models <- hr_models %>%
#  rbind(c("interactive_vs_controllability", as.numeric(hr_models[5, 2:16] / hr_models[1, 2:16])))


hr_models_long <- hr_models %>%
  #filter(model != "interactive_vs_controllability")%>%
  pivot_longer(
    !model,
    values_to = "BF",
    names_to = "timebin"
  )%>%
  mutate(
    time = as.numeric(timebin)-0.5,
    BF =as.numeric(BF)
  )

hr_models_long$model <- factor(hr_models_long$model, levels = c("controllability", "intensity", "interactive", "controllability_plus_interaction")) #"additive", "weighted", 

bayes_hr_plot <- ggplot(hr_models_long, aes(x = time, y= BF,  group = model, linetype = model, color = model)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE", color = NA)+
  geom_hline(yintercept=1, color = "black")+ 
  geom_line(linewidth = 0.8)+
  scale_linetype_manual(name="Model",  values=c("solid", "solid", "dotdash", "dotdash"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model"))+
  scale_color_manual(name="Model", values = c("#239BE7","#3434AD","#38d4a2", "#1c5c5a"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model")) +
  scale_y_continuous(name = "Bayes Factor", trans='log2')+  #breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000)
  scale_x_continuous(name = "Time",limits = c(0.2,15), breaks = seq(0,15,2))+
  theme_classic()


# SCL Data ---------------------------------------------------------------------

scl_long <- readRDS("EDA_unified.RData") %>%
  rename(cue= condition) %>%
  pivot_longer(
    cols = starts_with("scl_"),
    values_to = "scl",
    names_prefix = "scl_",
    names_to = "time"
  )%>%
  group_by(ID, time, cue)%>%
  summarize(
    value = mean(scl, na.rm = T)
  )%>%
  mutate(
    ID = as.factor(ID),
    time = as.numeric(time),
    predictor_exclusive_controllability = case_when(cue == 1 ~ -1,
                                                    cue == 2 ~ 1,
                                                    cue == 3 ~ -1,
                                                    cue == 4 ~ 1),
    predictor_exclusive_intensity = case_when(cue == 1 ~ 1,
                                              cue == 2 ~ 1,
                                              cue == 3 ~ -1,
                                              cue == 4 ~ -1),
    predictor_additive = case_when(cue == 1 ~ 0,
                                   cue == 2 ~ 2,
                                   cue == 3 ~ -2,
                                   cue == 4 ~ 0),
    predictor_weighted_add = case_when(cue == 1 ~ -1,
                                       cue == 2 ~ 3,
                                       cue == 3 ~ -3,
                                       cue == 4 ~ 1),
    predictor_interactive = case_when(cue == 1 ~ -3, #-1 #-2
                                      cue == 2 ~ 5,  #4  #4
                                      cue == 3 ~ -4, #-3 #-3
                                      cue == 4 ~ 2)) #0  #2

# Prepare Dataframe for Bayes Factors
scl_models <- data.frame(matrix(ncol = 11, nrow = 4))
names <- c("model", paste0(1:10))
colnames(scl_models) <- names
scl_models$model <- c("controllability", "intensity", "interactive", "controllability_plus_interaction" )#"additive", "weighted"

for (timebin in 1:10) {
  
  # Model 1: Controllability
  set.seed(30)
  BFfull_exclusive_controllability <- generalTestBF(value ~ predictor_exclusive_controllability + ID, whichRandom="ID", data=scl_long%>% filter(time == timebin), 
                                                    neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability <- BFfull_exclusive_controllability[1]/BFfull_exclusive_controllability[2]
  extractBF(BF_controllability)$bf
  scl_models[1,timebin+1] <- extractBF(BF_controllability)$bf
  
  # Model 2: Intensity
  set.seed(30)
  BFfull_exclusive_intensity <- generalTestBF(value ~ predictor_exclusive_intensity + ID, whichRandom="ID", data=scl_long%>% filter(time == timebin), 
                                              neverExclude="ID", iterations=100000, whichModels = "all")
  BF_intensity <- BFfull_exclusive_intensity[1]/BFfull_exclusive_intensity[2]
  scl_models[2,timebin+1] <- extractBF(BF_intensity)$bf
  
  # # Model 3: Controllability + Intensity
  # set.seed(30)
  # BFfull_additive <- generalTestBF(value ~ predictor_additive + ID, whichRandom="ID", data= scl_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_additive <- BFfull_additive[1]/BFfull_additive[2]
  # 
  # #BF_additive_vs_controllability <- BFfull_additive[1]/BFfull_exclusive_controllability[1]
  # scl_models[3,timebin+1] <- extractBF(BF_additive)$bf
  # 
  # # Model 3: Controllability + Intensity Weighted
  # set.seed(30)
  # BFfull_weighted <- generalTestBF(value ~ predictor_weighted_add + ID, whichRandom="ID", data=scl_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_weighted <- BFfull_weighted[1]/BFfull_weighted[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # scl_models[4,timebin+1] <- extractBF(BF_weighted)$bf
  # 
  # # Model 5: Controllability * Intensity Interaction
  # set.seed(30)
  # BFfull_interactive <- generalTestBF(value ~ predictor_interactive + ID, whichRandom="ID", data=scl_long%>% filter(time == timebin), 
  #                                     neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_interactive <- BFfull_interactive[1]/BFfull_interactive[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # scl_models[5,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 6: Interaction Effect 
  set.seed(30)
  BFfull_interactive <- generalTestBF(value ~ predictor_exclusive_controllability * predictor_exclusive_intensity + ID, whichRandom="ID", data=scl_long%>% filter(time == timebin), 
                                      neverExclude="ID", iterations=100000, whichModels = "all")
  BF_interactive <- BFfull_interactive[7]/BFfull_interactive[8]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  scl_models[3,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 7: Controllability plus Interaction Effect 
  set.seed(30)
  BFfull_controllability_interaction <- generalTestBF(value ~ predictor_exclusive_controllability + predictor_exclusive_controllability:predictor_exclusive_intensity + ID, whichRandom="ID", data=scl_long%>% filter(time == timebin), 
                                                      neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability_interaction <- BFfull_controllability_interaction[3]/BFfull_controllability_interaction[4]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  scl_models[4,timebin+1] <- extractBF(BF_controllability_interaction)$bf
}

#scl_models <- scl_models %>%
#  rbind(c("interactive_vs_intensity", as.numeric(scl_models[5, 2:11] / scl_models[2, 2:11])))

scl_models_long <- scl_models %>%
  #filter(model != "interactive_vs_intensity")%>%
  pivot_longer(
    !model,
    values_to = "BF",
    names_to = "timebin"
  )%>%
  mutate(
    time = as.numeric(timebin)-0.5,
    BF = as.numeric(BF)
  )

scl_models_long$model <- factor(scl_models_long$model, levels = c("controllability", "intensity", "interactive", "controllability_plus_interaction"))#"additive", "weighted"

bayes_scl_plot <- ggplot(scl_models_long%>% filter(model != "interactive_vs_intensity"), aes(x = time, y= BF,  group = model, linetype = model, color = model)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE", color = NA)+
  geom_hline(yintercept=1, color = "black")+ 
  geom_line(linewidth = 0.8)+
  scale_linetype_manual(name="Model",  values=c("solid", "solid", "dotdash", "dotdash"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model"))+
  scale_color_manual(name="Model", values = c("#239BE7","#3434AD","#38d4a2", "#1c5c5a"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model")) +
  scale_y_continuous(name = "Bayes Factor", trans='log2')+  #breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000)
  scale_x_continuous(name = "Time",limits = c(0.2,10), breaks = seq(0,10,2))+
  theme_classic()

# Pupil Data ----------------------------------------------------------------------
# ggplot(pupil_long%>%group_by(cue, time)%>% summarize(value = mean(value, na.rm = T)), aes(x= time, y = value, group = cue, color=cue))+
#   geom_line()

pupil_long <- readRDS("ET_pupil_df.RData") %>%
  pivot_longer(
    cols = starts_with("dilation_"),
    values_to = "pd",
    names_prefix = "dilation_",
    names_to = "time"
  )%>%
  group_by(ID, time, condition)%>%
  summarize(
    value = mean(pd, na.rm = T)
  )%>%
  mutate(
    ID = as.factor(ID),
    cue = condition,
    time = as.numeric(time),
    predictor_exclusive_controllability = case_when(cue == 1 ~ -1,
                                                    cue == 2 ~ 1,
                                                    cue == 3 ~ -1,
                                                    cue == 4 ~ 1),
    predictor_exclusive_intensity = case_when(cue == 1 ~ 1,
                                              cue == 2 ~ 1,
                                              cue == 3 ~ -1,
                                              cue == 4 ~ -1),
    predictor_additive = case_when(cue == 1 ~ 0,
                                   cue == 2 ~ 2,
                                   cue == 3 ~ -2,
                                   cue == 4 ~ 0),
    predictor_weighted_add = case_when(cue == 1 ~ -1,
                                       cue == 2 ~ 3,
                                       cue == 3 ~ -3,
                                       cue == 4 ~ 1),
    predictor_interactive = case_when(cue == 1 ~ -3, #-1 #-2
                                      cue == 2 ~ 5,  #4  #4
                                      cue == 3 ~ -4, #-3 #-3
                                      cue == 4 ~ 2)) #0  #2

# Prepare Dataframe for Bayes Factors
pupil_models <- data.frame(matrix(ncol = 21, nrow = 4))
names <- c("model", paste0(1:20))
colnames(pupil_models) <- names
pupil_models$model <- c("controllability", "intensity", "interactive", "controllability_plus_interaction") #"additive",  "weighted", 

for (timebin in 1:20) {
  
  # Model 1: Controllability
  set.seed(30)
  BFfull_exclusive_controllability <- generalTestBF(value ~ predictor_exclusive_controllability + ID, whichRandom="ID", data=pupil_long%>% filter(time == timebin), 
                                                    neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability <- BFfull_exclusive_controllability[1]/BFfull_exclusive_controllability[2]
  extractBF(BF_controllability)$bf
  pupil_models[1,timebin+1] <- extractBF(BF_controllability)$bf
  
  # Model 2: Intensity
  set.seed(30)
  BFfull_exclusive_intensity <- generalTestBF(value ~ predictor_exclusive_intensity + ID, whichRandom="ID", data=pupil_long%>% filter(time == timebin), 
                                              neverExclude="ID", iterations=100000, whichModels = "all")
  BF_intensity <- BFfull_exclusive_intensity[1]/BFfull_exclusive_intensity[2]
  pupil_models[2,timebin+1] <- extractBF(BF_intensity)$bf
  
  # # Model 3: Controllability + Intensity
  # set.seed(30)
  # BFfull_additive <- generalTestBF(value ~ predictor_additive + ID, whichRandom="ID", data= pupil_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_additive <- BFfull_additive[1]/BFfull_additive[2]
  # 
  # #BF_additive_vs_controllability <- BFfull_additive[1]/BFfull_exclusive_controllability[1]
  # pupil_models[3,timebin+1] <- extractBF(BF_additive)$bf
  # 
  # # Model 3: Controllability + Intensity Weighted
  # set.seed(30)
  # BFfull_weighted <- generalTestBF(value ~ predictor_weighted_add + ID, whichRandom="ID", data=pupil_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_weighted <- BFfull_weighted[1]/BFfull_weighted[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # pupil_models[4,timebin+1] <- extractBF(BF_weighted)$bf
  # 
  # # Model 5: Controllability * Intensity Interaction
  # set.seed(30)
  # BFfull_interactive <- generalTestBF(value ~ predictor_interactive + ID, whichRandom="ID", data=pupil_long%>% filter(time == timebin), 
  #                                     neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_interactive <- BFfull_interactive[1]/BFfull_interactive[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # pupil_models[3,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 6: Interaction Effect 
  set.seed(30)
  BFfull_interactive <- generalTestBF(value ~ predictor_exclusive_controllability * predictor_exclusive_intensity + ID, whichRandom="ID", data=pupil_long%>% filter(time == timebin), 
                                      neverExclude="ID", iterations=100000, whichModels = "all")
  BF_interactive <- BFfull_interactive[7]/BFfull_interactive[8]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  pupil_models[3,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 7: Controllability plus Interaction Effect 
  set.seed(30)
  BFfull_controllability_interaction <- generalTestBF(value ~ predictor_exclusive_controllability + predictor_exclusive_controllability:predictor_exclusive_intensity + ID, whichRandom="ID", data=pupil_long%>% filter(time == timebin), 
                                                      neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability_interaction <- BFfull_controllability_interaction[3]/BFfull_controllability_interaction[4]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  pupil_models[4,timebin+1] <- extractBF(BF_controllability_interaction)$bf
}

#pupil_models <- pupil_models %>%
#  rbind(c("interactive_vs_controllability", as.numeric(pupil_models[5, 2:21] / pupil_models[1, 2:21])))

pupil_models_long <- pupil_models %>%
  #filter(model != "interactive_vs_controllability")%>%
  pivot_longer(
    !model,
    values_to = "BF",
    names_to = "timebin"
  )%>%
  mutate(
    time = as.numeric(timebin)/2-0.25,
    BF = as.numeric(BF)
  )

pupil_models_long$model <- factor(pupil_models_long$model, levels = c("controllability", "intensity", "interactive", "controllability_plus_interaction")) #"additive", "weighted", 

bayes_pupil_plot <- ggplot(pupil_models_long, aes(x = time, y= BF,  group = model, linetype = model, color = model)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE", color = NA)+
  geom_hline(yintercept=1, color = "black")+ 
  geom_line(linewidth = 0.8)+
  scale_linetype_manual(name="Model",  values=c("solid", "solid", "dotdash", "dotdash"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model"))+
  scale_color_manual(name="Model", values = c("#239BE7","#3434AD","#38d4a2", "#1c5c5a"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model")) +
  scale_y_continuous(name = "Bayes Factor", trans='log2')+  #breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000)
  scale_x_continuous(name = "Time",limits = c(0.2,10), breaks = seq(0,10,2))+
  theme_classic()

# Central Bias -----------------------------------------------------------------

cb_long <- read.csv2(paste0(path, "Data/Tobii/CB_Eye_Data_Mean.csv")) %>%
  rename(ID = vp) %>%
  mutate(
    ID = as.factor(ID),
    cue = case_when(ActPass == "Passive" & Threat == "High" ~ 1,
                    ActPass == "Active" & Threat == "High" ~ 2,
                    ActPass == "Passive" & Threat == "Low" ~ 3,
                    ActPass == "Active" & Threat == "Low" ~ 4),
    time = bin,
    predictor_exclusive_controllability = case_when(cue == 1 ~ -1,
                                                    cue == 2 ~ 1,
                                                    cue == 3 ~ -1,
                                                    cue == 4 ~ 1),
    predictor_exclusive_intensity = case_when(cue == 1 ~ 1,
                                              cue == 2 ~ 1,
                                              cue == 3 ~ -1,
                                              cue == 4 ~ -1),
    predictor_additive = case_when(cue == 1 ~ 0,
                                   cue == 2 ~ 2,
                                   cue == 3 ~ -2,
                                   cue == 4 ~ 0),
    predictor_weighted_add = case_when(cue == 1 ~ -1,
                                       cue == 2 ~ 3,
                                       cue == 3 ~ -3,
                                       cue == 4 ~ 1),
    predictor_interactive = case_when(cue == 1 ~ -3, #-1 #-2
                                      cue == 2 ~ 5,  #4  #4
                                      cue == 3 ~ -4, #-3 #-3
                                      cue == 4 ~ 2)) #0  #2

# Prepare Dataframe for Bayes Factors
cb_models <- data.frame(matrix(ncol = 9, nrow = 4))
names <- c("model", paste0(1:8))
colnames(cb_models) <- names
cb_models$model <-  c("controllability", "intensity", "interactive", "controllability_plus_interaction") #"additive", "weighted", 

for (timebin in 1:8) {
  
  # Model 1: Controllability
  set.seed(30)
  BFfull_exclusive_controllability <- generalTestBF(value ~ predictor_exclusive_controllability + ID, whichRandom="ID", data=cb_long%>% filter(time == timebin), 
                                                    neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability <- BFfull_exclusive_controllability[1]/BFfull_exclusive_controllability[2]
  extractBF(BF_controllability)$bf
  cb_models[1,timebin+1] <- extractBF(BF_controllability)$bf
  
  # Model 2: Intensity
  set.seed(30)
  BFfull_exclusive_intensity <- generalTestBF(value ~ predictor_exclusive_intensity + ID, whichRandom="ID", data=cb_long%>% filter(time == timebin), 
                                              neverExclude="ID", iterations=100000, whichModels = "all")
  BF_intensity <- BFfull_exclusive_intensity[1]/BFfull_exclusive_intensity[2]
  cb_models[2,timebin+1] <- extractBF(BF_intensity)$bf
  
  # # Model 3: Controllability + Intensity
  # set.seed(30)
  # BFfull_additive <- generalTestBF(value ~ predictor_additive + ID, whichRandom="ID", data= cb_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_additive <- BFfull_additive[1]/BFfull_additive[2]
  # 
  # #BF_additive_vs_controllability <- BFfull_additive[1]/BFfull_exclusive_controllability[1]
  # cb_models[3,timebin+1] <- extractBF(BF_additive)$bf
  # 
  # # Model 3: Controllability + Intensity Weighted
  # set.seed(30)
  # BFfull_weighted <- generalTestBF(value ~ predictor_weighted_add + ID, whichRandom="ID", data=cb_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_weighted <- BFfull_weighted[1]/BFfull_weighted[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # cb_models[4,timebin+1] <- extractBF(BF_weighted)$bf
  # 
  # # Model 5: Controllability * Intensity Interaction
  # set.seed(30)
  # BFfull_interactive <- generalTestBF(value ~ predictor_interactive + ID, whichRandom="ID", data=cb_long%>% filter(time == timebin), 
  #                                     neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_interactive <- BFfull_interactive[1]/BFfull_interactive[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # cb_models[3,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 6: Interaction Effect 
  set.seed(30)
  BFfull_interactive <- generalTestBF(value ~ predictor_exclusive_controllability * predictor_exclusive_intensity + ID, whichRandom="ID", data=cb_long%>% filter(time == timebin), 
                                      neverExclude="ID", iterations=100000, whichModels = "all")
  BF_interactive <- BFfull_interactive[7]/BFfull_interactive[8]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  cb_models[3,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 7: Controllability plus Interaction Effect 
  set.seed(30)
  BFfull_controllability_interaction <- generalTestBF(value ~ predictor_exclusive_controllability + predictor_exclusive_controllability:predictor_exclusive_intensity + ID, whichRandom="ID", data=cb_long%>% filter(time == timebin), 
                                                      neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability_interaction <- BFfull_controllability_interaction[3]/BFfull_controllability_interaction[4]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  cb_models[4,timebin+1] <- extractBF(BF_controllability_interaction)$bf
}

#cb_models <- cb_models %>%
#  rbind(c("interactive_vs_controllability", as.numeric(cb_models[5, 2:9] / cb_models[1, 2:9])))

cb_models_long <- cb_models %>%
  #filter(model!= "interactive_vs_controllability")%>%
  pivot_longer(
    !model,
    values_to = "BF",
    names_to = "timebin"
  )%>%
  mutate(
    time = as.numeric(timebin)+1.5,
    BF = as.numeric(BF)
  )

cb_models_long$model <- factor(cb_models_long$model, levels = c("controllability", "intensity", "interactive", "controllability_plus_interaction")) #"additive", "weighted",

bayes_cb_plot <- ggplot(cb_models_long, aes(x = time, y= BF,  group = model, linetype = model, color = model)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE", color = NA)+
  geom_hline(yintercept=1, color = "black")+ 
  geom_line(linewidth = 0.8)+
  scale_linetype_manual(name="Model",  values=c("solid", "solid", "dotdash", "dotdash"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model"))+
  scale_color_manual(name="Model", values = c("#239BE7","#3434AD","#38d4a2", "#1c5c5a"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model")) +
  scale_y_continuous(name = "Bayes Factor", trans='log2')+  #breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000)
  scale_x_continuous(name = "Time",limits = c(0.2,10), breaks = seq(0,10,2))+
  theme_classic()

# Fix Dur -----------------------------------------------------------------

fd_long <- read.csv2(paste0(path, "Data/Tobii/FixDur_Eye_Data_Mean.csv")) %>%
  rename(ID = vp) %>%
  mutate(
    ID = as.factor(ID),
    cue = case_when(ActPass == "Passive" & Threat == "High" ~ 1,
                    ActPass == "Active" & Threat == "High" ~ 2,
                    ActPass == "Passive" & Threat == "Low" ~ 3,
                    ActPass == "Active" & Threat == "Low" ~ 4),
    time = bin,
    predictor_exclusive_controllability = case_when(cue == 1 ~ -1,
                                                    cue == 2 ~ 1,
                                                    cue == 3 ~ -1,
                                                    cue == 4 ~ 1),
    predictor_exclusive_intensity = case_when(cue == 1 ~ 1,
                                              cue == 2 ~ 1,
                                              cue == 3 ~ -1,
                                              cue == 4 ~ -1),
    predictor_additive = case_when(cue == 1 ~ 0,
                                   cue == 2 ~ 2,
                                   cue == 3 ~ -2,
                                   cue == 4 ~ 0),
    predictor_weighted_add = case_when(cue == 1 ~ -1,
                                       cue == 2 ~ 3,
                                       cue == 3 ~ -3,
                                       cue == 4 ~ 1),
    predictor_interactive = case_when(cue == 1 ~ -3, #-1 #-2
                                      cue == 2 ~ 5,  #4  #4
                                      cue == 3 ~ -4, #-3 #-3
                                      cue == 4 ~ 2)) #0  #2

fd_long %>% 
  group_by(cue)%>%
  summarize(fd = mean(value),
            sd = sd(value))

# Prepare Dataframe for Bayes Factors
fd_models <- data.frame(matrix(ncol = 9, nrow = 4))
names <- c("model", paste0(1:8))
colnames(fd_models) <- names
fd_models$model <- c("controllability", "intensity", "interactive", "controllability_plus_interaction")#"additive", "weighted",

for (timebin in 1:8) {
  
  # Model 1: Controllability
  set.seed(30)
  BFfull_exclusive_controllability <- generalTestBF(value ~ predictor_exclusive_controllability + ID, whichRandom="ID", data=fd_long%>% filter(time == timebin), 
                                                    neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability <- BFfull_exclusive_controllability[1]/BFfull_exclusive_controllability[2]
  extractBF(BF_controllability)$bf
  fd_models[1,timebin+1] <- extractBF(BF_controllability)$bf
  
  # Model 2: Intensity
  set.seed(30)
  BFfull_exclusive_intensity <- generalTestBF(value ~ predictor_exclusive_intensity + ID, whichRandom="ID", data=fd_long%>% filter(time == timebin), 
                                              neverExclude="ID", iterations=100000, whichModels = "all")
  BF_intensity <- BFfull_exclusive_intensity[1]/BFfull_exclusive_intensity[2]
  fd_models[2,timebin+1] <- extractBF(BF_intensity)$bf
  
  # # Model 3: Controllability + Intensity
  # set.seed(30)
  # BFfull_additive <- generalTestBF(value ~ predictor_additive + ID, whichRandom="ID", data= fd_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_additive <- BFfull_additive[1]/BFfull_additive[2]
  # 
  # #BF_additive_vs_controllability <- BFfull_additive[1]/BFfull_exclusive_controllability[1]
  # fd_models[3,timebin+1] <- extractBF(BF_additive)$bf
  
  # # Model 3: Controllability x Intensity
  # set.seed(30)
  # BFfull_weighted <- generalTestBF(value ~ predictor_weighted_add + ID, whichRandom="ID", data=fd_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_weighted <- BFfull_weighted[1]/BFfull_weighted[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # fd_models[4,timebin+1] <- extractBF(BF_weighted)$bf
  # 
  # # Model 5: Controllability * Intensity Interaction
  # set.seed(30)
  # BFfull_interactive <- generalTestBF(value ~ predictor_interactive + ID, whichRandom="ID", data=fd_long%>% filter(time == timebin), 
  #                                     neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_interactive <- BFfull_interactive[1]/BFfull_interactive[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # fd_models[3,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 6: Interaction Effect 
  set.seed(30)
  BFfull_interactive <- generalTestBF(value ~ predictor_exclusive_controllability * predictor_exclusive_intensity + ID, whichRandom="ID", data=fd_long%>% filter(time == timebin), 
                                      neverExclude="ID", iterations=100000, whichModels = "all")
  BF_interactive <- BFfull_interactive[7]/BFfull_interactive[8]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  fd_models[3,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 7: Controllability plus Interaction Effect 
  set.seed(30)
  BFfull_controllability_interaction <- generalTestBF(value ~ predictor_exclusive_controllability + predictor_exclusive_controllability:predictor_exclusive_intensity + ID, whichRandom="ID", data=fd_long%>% filter(time == timebin), 
                                                      neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability_interaction <- BFfull_controllability_interaction[3]/BFfull_controllability_interaction[4]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  fd_models[4,timebin+1] <- extractBF(BF_controllability_interaction)$bf
}

#fd_models <- fd_models %>%
#  rbind(c("interactive_vs_controllability", as.numeric(fd_models[5, 2:9] / fd_models[1, 2:9])))

fd_models_long <- fd_models %>%
  #filter(model != "interactive_vs_controllability")%>%
  pivot_longer(
    !model,
    values_to = "BF",
    names_to = "timebin"
  )%>%
  mutate(
    time = as.numeric(timebin)+1.5,
    BF = as.numeric(BF)
  )

fd_models_long$model <- factor(fd_models_long$model, levels = c("controllability", "intensity", "interactive", "controllability_plus_interaction"))

bayes_fd_plot <- ggplot(fd_models_long, aes(x = time, y= BF,  group = model, linetype = model, color = model)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE", color = NA)+
  geom_hline(yintercept=1, color = "black")+ 
  geom_line(linewidth = 0.8)+
  scale_linetype_manual(name="Model",  values=c("solid", "solid", "dotdash", "dotdash"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model"))+
  scale_color_manual(name="Model", values = c("#239BE7","#3434AD","#38d4a2", "#1c5c5a"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model")) +
  scale_y_continuous(name = "Bayes Factor", trans='log2')+  #breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000)
  scale_x_continuous(name = "Time",limits = c(0.2,10), breaks = seq(0,10,2))+
  theme_classic()


# Fix Num -----------------------------------------------------------------

fn_long <- read.csv2(paste0(path, "Data/Tobii/FixNum_Eye_Data_Mean.csv")) %>%
  rename(ID = vp) %>%
  mutate(
    ID = as.factor(ID),
    cue = case_when(ActPass == "Passive" & Threat == "High" ~ 1,
                    ActPass == "Active" & Threat == "High" ~ 2,
                    ActPass == "Passive" & Threat == "Low" ~ 3,
                    ActPass == "Active" & Threat == "Low" ~ 4),
    time = bin,
    predictor_exclusive_controllability = case_when(cue == 1 ~ -1,
                                                    cue == 2 ~ 1,
                                                    cue == 3 ~ -1,
                                                    cue == 4 ~ 1),
    predictor_exclusive_intensity = case_when(cue == 1 ~ 1,
                                              cue == 2 ~ 1,
                                              cue == 3 ~ -1,
                                              cue == 4 ~ -1),
    predictor_additive = case_when(cue == 1 ~ 0,
                                   cue == 2 ~ 2,
                                   cue == 3 ~ -2,
                                   cue == 4 ~ 0),
    predictor_weighted_add = case_when(cue == 1 ~ -1,
                                       cue == 2 ~ 3,
                                       cue == 3 ~ -3,
                                       cue == 4 ~ 1),
    predictor_interactive = case_when(cue == 1 ~ -3, #-1 #-2
                                      cue == 2 ~ 5,  #4  #4
                                      cue == 3 ~ -4, #-3 #-3
                                      cue == 4 ~ 2)) #0  #2

fn_long %>% 
  group_by(cue)%>%
  summarize(fn = mean(value),
            sd = sd(value))

# Prepare Dataframe for Bayes Factors
fn_models <- data.frame(matrix(ncol = 9, nrow = 4))
names <- c("model", paste0(1:8))
colnames(fn_models) <- names
fn_models$model <- c("controllability", "intensity", "interactive", "controllability_plus_interaction")

for (timebin in 1:8) {
  
  # Model 1: Controllability
  set.seed(30)
  BFfull_exclusive_controllability <- generalTestBF(value ~ predictor_exclusive_controllability + ID, whichRandom="ID", data=fn_long%>% filter(time == timebin), 
                                                    neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability <- BFfull_exclusive_controllability[1]/BFfull_exclusive_controllability[2]
  extractBF(BF_controllability)$bf
  fn_models[1,timebin+1] <- extractBF(BF_controllability)$bf
  
  # Model 2: Intensity
  set.seed(30)
  BFfull_exclusive_intensity <- generalTestBF(value ~ predictor_exclusive_intensity + ID, whichRandom="ID", data=fn_long%>% filter(time == timebin), 
                                              neverExclude="ID", iterations=100000, whichModels = "all")
  BF_intensity <- BFfull_exclusive_intensity[1]/BFfull_exclusive_intensity[2]
  fn_models[2,timebin+1] <- extractBF(BF_intensity)$bf
  
  # # Model 3: Controllability + Intensity
  # set.seed(30)
  # BFfull_additive <- generalTestBF(value ~ predictor_additive + ID, whichRandom="ID", data= fn_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_additive <- BFfull_additive[1]/BFfull_additive[2]
  # 
  # #BF_additive_vs_controllability <- BFfull_additive[1]/BFfull_exclusive_controllability[1]
  # fn_models[3,timebin+1] <- extractBF(BF_additive)$bf
  # 
  # # Model 3: Controllability x Intensity
  # set.seed(30)
  # BFfull_weighted <- generalTestBF(value ~ predictor_weighted_add + ID, whichRandom="ID", data=fn_long%>% filter(time == timebin), 
  #                                  neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_weighted <- BFfull_weighted[1]/BFfull_weighted[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # fn_models[4,timebin+1] <- extractBF(BF_weighted)$bf
  # 
  # # Model 5: Controllability * Intensity Interaction
  # set.seed(30)
  # BFfull_interactive <- generalTestBF(value ~ predictor_interactive + ID, whichRandom="ID", data=fn_long%>% filter(time == timebin), 
  #                                     neverExclude="ID", iterations=100000, whichModels = "all")
  # BF_interactive <- BFfull_interactive[1]/BFfull_interactive[2]
  # 
  # #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  # fn_models[3,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 6: Interaction Effect 
  set.seed(30)
  BFfull_interactive <- generalTestBF(value ~ predictor_exclusive_controllability * predictor_exclusive_intensity + ID, whichRandom="ID", data=fn_long%>% filter(time == timebin), 
                                      neverExclude="ID", iterations=100000, whichModels = "all")
  BF_interactive <- BFfull_interactive[7]/BFfull_interactive[8]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  fn_models[3,timebin+1] <- extractBF(BF_interactive)$bf
  
  # Model 7: Controllability plus Interaction Effect 
  set.seed(30)
  BFfull_controllability_interaction <- generalTestBF(value ~ predictor_exclusive_controllability + predictor_exclusive_controllability:predictor_exclusive_intensity + ID, whichRandom="ID", data=fn_long%>% filter(time == timebin), 
                                                      neverExclude="ID", iterations=100000, whichModels = "all")
  BF_controllability_interaction <- BFfull_controllability_interaction[3]/BFfull_controllability_interaction[4]
  
  #BF_interactive_vs_controllability <- BFfull_weighted[1]/BFfull_exclusive_controllability[1] 
  fn_models[4,timebin+1] <- extractBF(BF_controllability_interaction)$bf
}

#fn_models <- fn_models %>%
#  rbind(c("interactive_vs_controllability", as.numeric(fn_models[5, 2:9] / fn_models[1, 2:9])))

fn_models_long <- fn_models %>%
  #filter(model != "interactive_vs_controllability")%>%
  pivot_longer(
    !model,
    values_to = "BF",
    names_to = "timebin"
  )%>%
  mutate(
    time = as.numeric(timebin)+1.5,
    BF = as.numeric(BF)
  )

fn_models_long$model <- factor(fn_models_long$model, levels = c("controllability", "intensity", "interactive", "controllability_plus_interaction"))

bayes_fn_plot <- ggplot(fn_models_long, aes(x = time, y= BF,  group = model, linetype = model, color = model)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE", color = NA)+
  geom_hline(yintercept=1, color = "black")+ 
  geom_line(linewidth = 0.8)+
  scale_linetype_manual(name="Model",  values=c("solid", "solid", "dotdash", "dotdash"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model"))+
  scale_color_manual(name="Model", values = c("#239BE7","#3434AD","#38d4a2", "#1c5c5a"),labels = c("Controllability Main Effect Model","Intensity Main Effect Model","Full Interactive Model", "Controllability plus Interaction Model")) +
  scale_y_continuous(name = "Bayes Factor", trans='log2')+  #breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000)
  scale_x_continuous(name = "Time",limits = c(0.2,10), breaks = seq(0,10,2))+
  theme_classic()



# Plot Contrasts ---------------------------------------------------------------

contrasts_controllability_plot <- ggplot(alpha_long%>% group_by(cue) %>% summarize(predictor = mean(predictor_exclusive_controllability))%>% mutate(cue = case_when(cue == 1 ~ "HI passive",
                                                                                                                                                                    cue == 2 ~ "HI active",
                                                                                                                                                                    cue == 3 ~ "LI passive",
                                                                                                                                                                    cue == 4 ~ "LI active")), aes(x = cue, y = predictor)) +
  geom_bar(stat="identity")+
  #scale_y_continuous(name = "Bayes Factor", trans='log2', breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000))+
  scale_x_discrete(name = "Condition")+
  geom_hline(yintercept=0, color = "black")+ 
  ggtitle("Main Effect Model Controllability") +
  theme_classic()

contrasts_intensity_plot <- ggplot(alpha_long%>% group_by(cue) %>% summarize(predictor = mean(predictor_exclusive_intensity))%>% mutate(cue = case_when(cue == 1 ~ "HI passive",
                                                                                                                                                        cue == 2 ~ "HI active",
                                                                                                                                                        cue == 3 ~ "LI passive",
                                                                                                                                                        cue == 4 ~ "LI active")), aes(x = cue, y = predictor)) +
  geom_bar(stat="identity")+
  #scale_y_continuous(name = "Bayes Factor", trans='log2', breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000))+
  scale_x_discrete(name = "Condition")+
  geom_hline(yintercept=0, color = "black")+ 
  ggtitle("Main Effect Model Intensity") +
  theme_classic()

# contrasts_additive_plot <- ggplot(alpha_long%>% group_by(cue) %>% summarize(predictor = mean(predictor_additive))%>% mutate(cue = case_when(cue == 1 ~ "HI passive",
#                                                                                                                                             cue == 2 ~ "HI active",
#                                                                                                                                             cue == 3 ~ "LI passive",
#                                                                                                                                             cue == 4 ~ "LI active")), aes(x = cue, y = predictor)) +
#   geom_bar(stat="identity")+
#   #scale_y_continuous(name = "Bayes Factor", trans='log2', breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000))+
#   scale_x_discrete(name = "Condition")+
#   geom_hline(yintercept=0, color = "black")+ 
#   ggtitle("Additive Model Controllability + Intensity") +
#   theme_classic()
# 
# contrasts_additive_weighted_plot <- ggplot(alpha_long%>% group_by(cue) %>% summarize(predictor = mean(predictor_weighted_add))%>% mutate(cue = case_when(cue == 1 ~ "HI passive",
#                                                                                                                                                          cue == 2 ~ "HI active",
#                                                                                                                                                          cue == 3 ~ "LI passive",
#                                                                                                                                                          cue == 4 ~ "LI active")), aes(x = cue, y = predictor)) +
#   geom_bar(stat="identity")+
#   #scale_y_continuous(name = "Bayes Factor", trans='log2', breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000))+
#   scale_x_discrete(name = "Condition")+
#   geom_hline(yintercept=0, color = "black")+ 
#   ggtitle("Weighted Additive Model 2*Controllability + 1*Intensity") +
#   theme_classic()
# 
# contrasts_interactive_plot <- ggplot(alpha_long%>% group_by(cue) %>% summarize(predictor = mean(predictor_interactive))%>% mutate(cue = case_when(cue == 1 ~ "HI passive",
#                                                                                                                                                   cue == 2 ~ "HI active",
#                                                                                                                                                   cue == 3 ~ "LI passive",
#                                                                                                                                                   cue == 4 ~ "LI active")), aes(x = cue, y = predictor)) +
#   geom_bar(stat="identity")+
#   #scale_y_continuous(name = "Bayes Factor", trans='log2', breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000,10000))+
#   scale_x_discrete(name = "Condition")+
#   geom_hline(yintercept=0, color = "black")+ 
#   ggtitle("Interactive Model") +
#   theme_classic()

#contrasts_plots <- contrasts_controllability_plot + contrasts_intensity_plot + contrasts_additive_plot + contrasts_additive_weighted_plot + contrasts_interactive_plot

