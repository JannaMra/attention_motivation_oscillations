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

questionnaires = suppressMessages( #avoid name repair message (affected columns will be deselected anyway)
  readxl::read_excel("Questionnaires.xlsx" %>% paste0(path, "Data/", .),
                     sheet = "All_Questionnaires_Overview", col_names=T))

questionnaires %>% summarise(
  mean_age = mean(Age, na.rm = T),
  sd_age = sd(Age, na.rm = T),
  mean_pain_mA_block1 = mean(PR01_01, na.rm = T),
  sd_pain_mA_block1 = sd(PR01_01, na.rm = T),
  mean_pain_mA_block2 = mean(PR02_01, na.rm = T),
  sd_pain_mA_block2 = sd(PR02_01, na.rm = T),
  mean_pain_li_pre = mean(PR03_01, na.rm = T),
  sd_pain_li_pre = sd(PR03_01, na.rm = T),
  mean_pain_li_post = mean(PR04_01, na.rm = T),
  sd_pain_li_post = sd(PR04_01, na.rm = T),
)

questionnaires %>% 
  group_by(Gender) %>%
  summarise(
    n = n()
)

questionnaires %>% 
  group_by(Handedness) %>%
  summarise(
    n = n()
  )

