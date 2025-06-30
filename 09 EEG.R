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

library(tidyverse)
library(ez)
library(apa)
library(readxl)


# paths
path <- paste0(getwd(),"/")

# Exclusions:
exclusions <- c("vp05","vp22","vp30","vp51")
eye.invalid.bl <- c("vp04", "vp08", "vp09", "vp20", "vp31", "vp46")

#ANOVA 500-1000 ms after cue onset
alpha1 <- read_excel("Data/EEG/bins1/TFO_alphas_Pz7_2.xls")
alpha1_long <- alpha1 %>%
  pivot_longer(High_Passive:Low_Active, 
               names_to = "ActPass", values_to = "value")%>%
  mutate(Intensity = ifelse(str_detect(ActPass,"High"),"High","Low"))%>%
  mutate(ActPass = gsub("High_","",ActPass))%>%
  mutate(ActPass = gsub("Low_","",ActPass))
alpha1_anova = ezANOVA(data=alpha1_long, dv=.(value), wid=.(ID), within=.(ActPass,Intensity), detailed=T, type=3)
alpha1_anova%>%
  apa::anova_apa(force_sph_corr=T)

print(alpha1_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # ActPass
print(alpha1_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # Intensity
print(alpha1_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # Interaction

#ANOVA 2500-3000 ms after cue onset
alpha2 <- read_excel("Data/EEG/bins1/TFO_alphas_Pz7_6.xls")
alpha2_long <- alpha2 %>%
  pivot_longer(High_Passive:Low_Active, 
               names_to = "ActPass", values_to = "value")%>%
  mutate(Intensity = ifelse(str_detect(ActPass,"High"),"High","Low"))%>%
  mutate(ActPass = gsub("High_","",ActPass))%>%
  mutate(ActPass = gsub("Low_","",ActPass))
alpha2_anova = ezANOVA(data=alpha2_long, dv=.(value), wid=.(ID), within=.(ActPass,Intensity), detailed=T, type=3)

alpha2_anova%>%
  apa::anova_apa(force_sph_corr=T)

print(alpha2_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # ActPass
print(alpha2_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # Intensity
print(alpha2_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # Interaction

#ANOVA 9500-10000 ms after cue onset
alpha3 <- read_excel("Data/EEG/bins1/TFO_alphas_Pz7_20.xls")
alpha3_long <- alpha3 %>%
  pivot_longer(High_Passive:Low_Active, 
               names_to = "ActPass", values_to = "value")%>%
  mutate(Intensity = ifelse(str_detect(ActPass,"High"),"High","Low"))%>%
  mutate(ActPass = gsub("High_","",ActPass))%>%
  mutate(ActPass = gsub("Low_","",ActPass))
alpha3_anova = ezANOVA(data=alpha3_long, dv=.(value), wid=.(ID), within=.(ActPass,Intensity), detailed=T, type=3)%>%
  apa::anova_apa(force_sph_corr=T)

#ANOVA across all bins 

alpha.df <- data.frame()
for (bin in 1:20) {
  alpha <- read_excel(paste0("Data/EEG/bins1/TFO_alphas_Pz7_",bin,".xls"))
  alpha <- alpha %>% mutate(time = as.numeric(bin))
  alpha.df <- rbind(alpha.df,alpha)
}

alpha.df_long <- alpha.df %>%
  pivot_longer(High_Passive:Low_Active, 
               names_to = "ActPass", values_to = "value")%>%
  mutate(Intensity = ifelse(str_detect(ActPass,"High"),"High","Low"))%>%
  mutate(ActPass = gsub("High_","",ActPass))%>%
  mutate(ActPass = gsub("Low_","",ActPass)) %>%
  mutate(time = as.factor(time), ActPass = as.factor(ActPass), Intensity = as.factor(Intensity))
write.csv2(alpha.df_long,"alpha_long.csv",row.names=FALSE,quote=FALSE)

read.csv2("alpha_long.csv")

alpha_anova = ezANOVA(data=alpha.df_long, dv=.(value), wid=.(ID), within=.(ActPass,Intensity,time), detailed=T, type=3)

alpha_anova %>%
  apa::anova_apa(force_sph_corr=T)

print(alpha_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # ActPass
print(alpha_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # Intensity
print(alpha_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # time
print(alpha_anova$ANOVA[5,] %>% partial_eta_squared_ci()) # ActPass * Intensity
print(alpha_anova$ANOVA[6,] %>% partial_eta_squared_ci()) # ActPass * Time
print(alpha_anova$ANOVA[7,] %>% partial_eta_squared_ci()) # Intensity * Time
print(alpha_anova$ANOVA[8,] %>% partial_eta_squared_ci()) # Time * ActPass * Threat


########################################################################################
### High Alpha #########################################################################
########################################################################################

#ANOVA 500-1000 ms after cue onset
alpha1_high <- read_excel("Data/EEG/bins1/TFO_high_alphas_Pz7_2.xls")
alpha1_high_long <- alpha1_high %>%
  pivot_longer(High_Passive:Low_Active, 
               names_to = "ActPass", values_to = "value")%>%
  mutate(Intensity = ifelse(str_detect(ActPass,"High"),"High","Low"))%>%
  mutate(ActPass = gsub("High_","",ActPass))%>%
  mutate(ActPass = gsub("Low_","",ActPass))
alpha1_high_anova = ezANOVA(data=alpha1_high_long, dv=.(value), wid=.(ID), within=.(ActPass,Intensity), detailed=T, type=3)
alpha1_high_anova%>%
  apa::anova_apa(force_sph_corr=T)

print(alpha1_high_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # ActPass
print(alpha1_high_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # Intensity
print(alpha1_high_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # Interaction

#ANOVA 2500-3000 ms after cue onset
alpha2_high <- read_excel("Data/EEG/bins1/TFO_high_alphas_Pz7_6.xls")
alpha2_high_long <- alpha2_high %>%
  pivot_longer(High_Passive:Low_Active, 
               names_to = "ActPass", values_to = "value")%>%
  mutate(Intensity = ifelse(str_detect(ActPass,"High"),"High","Low"))%>%
  mutate(ActPass = gsub("High_","",ActPass))%>%
  mutate(ActPass = gsub("Low_","",ActPass))
alpha2_high_anova = ezANOVA(data=alpha2_high_long, dv=.(value), wid=.(ID), within=.(ActPass,Intensity), detailed=T, type=3)

alpha2_high_anova%>%
  apa::anova_apa(force_sph_corr=T)

print(alpha2_high_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # ActPass
print(alpha2_high_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # Intensity
print(alpha2_high_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # Interaction

#ANOVA 9500-10000 ms after cue onset
alpha3_high <- read_excel("Data/EEG/bins1/TFO_high_alphas_Pz7_20.xls")
alpha3_high_long <- alpha3_high %>%
  pivot_longer(High_Passive:Low_Active, 
               names_to = "ActPass", values_to = "value")%>%
  mutate(Intensity = ifelse(str_detect(ActPass,"High"),"High","Low"))%>%
  mutate(ActPass = gsub("High_","",ActPass))%>%
  mutate(ActPass = gsub("Low_","",ActPass))
alpha3_high_anova = ezANOVA(data=alpha3_high_long, dv=.(value), wid=.(ID), within=.(ActPass,Intensity), detailed=T, type=3)

alpha3_high_anova%>%
  apa::anova_apa(force_sph_corr=T)

print(alpha3_high_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # ActPass
print(alpha3_high_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # Intensity
print(alpha3_high_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # Interaction

#ANOVA across all bins high alpha

high.alpha.df <- data.frame()
for (bin in 1:20) {
  alpha <- read_excel(paste0("Data/EEG/bins1/TFO_high_alphas_Pz7_",bin,".xls"))
  alpha <- alpha %>% mutate(time = as.numeric(bin))
  high.alpha.df <- rbind(high.alpha.df,alpha)
}

high.alpha.df_long <- high.alpha.df %>%
  pivot_longer(High_Passive:Low_Active, 
               names_to = "ActPass", values_to = "value")%>%
  mutate(Intensity = ifelse(str_detect(ActPass,"High"),"High","Low"))%>%
  mutate(ActPass = gsub("High_","",ActPass))%>%
  mutate(ActPass = gsub("Low_","",ActPass)) %>%
  mutate(time = as.factor(time), ActPass = as.factor(ActPass), Intensity = as.factor(Intensity))
write.csv2(high.alpha.df_long,"high_alpha_long.csv",row.names=FALSE,quote=FALSE)

high_alpha_anova = ezANOVA(data=high.alpha.df_long, dv=.(value), wid=.(ID), within=.(ActPass,Intensity,time), detailed=T, type=3)

high_alpha_anova %>%
  apa::anova_apa(force_sph_corr=T)

print(high_alpha_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # ActPass
print(high_alpha_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # Intensity
print(high_alpha_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # time
print(high_alpha_anova$ANOVA[5,] %>% partial_eta_squared_ci()) # ActPass * Intensity
print(high_alpha_anova$ANOVA[6,] %>% partial_eta_squared_ci()) # ActPass * Time
print(high_alpha_anova$ANOVA[7,] %>% partial_eta_squared_ci()) # Intensity * Time
print(high_alpha_anova$ANOVA[8,] %>% partial_eta_squared_ci()) # Time * ActPass * Threat

########################################################################################
### Low Alpha #########################################################################
########################################################################################

#ANOVA 500-1000 ms after cue onset
alpha1_low <- read_excel("Data/EEG/bins1/TFO_low_alphas_Pz7_2.xls")
alpha1_low_long <- alpha1_low %>%
  pivot_longer(High_Passive:Low_Active, 
               names_to = "ActPass", values_to = "value")%>%
  mutate(Intensity = ifelse(str_detect(ActPass,"High"),"High","Low"))%>%
  mutate(ActPass = gsub("High_","",ActPass))%>%
  mutate(ActPass = gsub("Low_","",ActPass))
alpha1_low_anova = ezANOVA(data=alpha1_low_long, dv=.(value), wid=.(ID), within=.(ActPass,Intensity), detailed=T, type=3)
alpha1_low_anova%>%
  apa::anova_apa(force_sph_corr=T)

print(alpha1_low_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # ActPass
print(alpha1_low_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # Intensity
print(alpha1_low_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # Interaction

#ANOVA 2500-3000 ms after cue onset
alpha2_low <- read_excel("Data/EEG/bins1/TFO_low_alphas_Pz7_6.xls")
alpha2_low_long <- alpha2_low %>%
  pivot_longer(High_Passive:Low_Active, 
               names_to = "ActPass", values_to = "value")%>%
  mutate(Intensity = ifelse(str_detect(ActPass,"High"),"High","Low"))%>%
  mutate(ActPass = gsub("High_","",ActPass))%>%
  mutate(ActPass = gsub("Low_","",ActPass))
alpha2_low_anova = ezANOVA(data=alpha2_low_long, dv=.(value), wid=.(ID), within=.(ActPass,Intensity), detailed=T, type=3)

alpha2_low_anova%>%
  apa::anova_apa(force_sph_corr=T)

print(alpha2_low_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # ActPass
print(alpha2_low_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # Intensity
print(alpha2_low_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # Interaction

#ANOVA 9500-10000 ms after cue onset
alpha3_low <- read_excel("Data/EEG/bins1/TFO_low_alphas_Pz7_20.xls")
alpha3_low_long <- alpha3_low %>%
  pivot_longer(High_Passive:Low_Active, 
               names_to = "ActPass", values_to = "value")%>%
  mutate(Intensity = ifelse(str_detect(ActPass,"High"),"High","Low"))%>%
  mutate(ActPass = gsub("High_","",ActPass))%>%
  mutate(ActPass = gsub("Low_","",ActPass))
alpha3_low_anova = ezANOVA(data=alpha3_low_long, dv=.(value), wid=.(ID), within=.(ActPass,Intensity), detailed=T, type=3)

alpha3_low_anova%>%
  apa::anova_apa(force_sph_corr=T)

print(alpha3_low_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # ActPass
print(alpha3_low_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # Intensity
print(alpha3_low_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # Interaction

#ANOVA across all bins low alpha

low.alpha.df <- data.frame()
for (bin in 1:20) {
  alpha <- read_excel(paste0("Data/EEG/bins1/TFO_low_alphas_Pz7_",bin,".xls"))
  alpha <- alpha %>% mutate(time = as.numeric(bin))
  low.alpha.df <- rbind(low.alpha.df,alpha)
}

low.alpha.df_long <- low.alpha.df %>%
  pivot_longer(High_Passive:Low_Active, 
               names_to = "ActPass", values_to = "value")%>%
  mutate(Intensity = ifelse(str_detect(ActPass,"High"),"High","Low"))%>%
  mutate(ActPass = gsub("High_","",ActPass))%>%
  mutate(ActPass = gsub("Low_","",ActPass)) %>%
  mutate(time = as.factor(time), ActPass = as.factor(ActPass), Intensity = as.factor(Intensity))
write.csv2(low.alpha.df_long,"low_alpha_long.csv",row.names=FALSE,quote=FALSE)

low_alpha_anova = ezANOVA(data=low.alpha.df_long, dv=.(value), wid=.(ID), within=.(ActPass,Intensity,time), detailed=T, type=3)

low_alpha_anova %>%
  apa::anova_apa(force_sph_corr=T)

print(low_alpha_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # ActPass
print(low_alpha_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # Intensity
print(low_alpha_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # time
print(low_alpha_anova$ANOVA[5,] %>% partial_eta_squared_ci()) # ActPass * Intensity
print(low_alpha_anova$ANOVA[6,] %>% partial_eta_squared_ci()) # ActPass * Time
print(low_alpha_anova$ANOVA[7,] %>% partial_eta_squared_ci()) # Intensity * Time
print(low_alpha_anova$ANOVA[8,] %>% partial_eta_squared_ci()) # Time * ActPass * Threat

