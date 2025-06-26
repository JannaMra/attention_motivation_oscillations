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
#install.packages("apa")
library(tidyverse)
library(ez)
library(apa)

path <- paste0(getwd(),"/")

# Onsets laden
msg <- read.table(paste(path,"Data/Tobii/Messages.txt",sep=""))
names(msg) <- c("vp","trial","time")
fixa <- read.table(paste(path,"Data/Tobii/Fixations.txt", sep=""))

# Determine which subjects should be analyzed
#vpn <- c(paste("tfo",ifelse(1:16<10,"0",""),1:16,sep=""))
vpn = fixa$subject %>% unique() %>% sort() %>% as.character() #all subjects in fixations
vpn.n = length(vpn)
vpn <-  vpn[!(vpn %in% exclusions)]
vpn <-  vpn[!(vpn %in% eye.invalid.bl)]
vpn.n = length(vpn)


# Merge
All_Eye_Data <- data.frame()

for (vp in vpn) { 
  
  #vp <- "vp01"
  vpcode <- vp
  print(vpcode)
  
  # read in data
  All_Eye_Data_VP <- read.csv2(paste(path.prot,vpcode,"_et.csv",sep=""), header = TRUE, sep=";", dec=",", check.names=F)
  All_Eye_Data_VP$vp <- vpcode
  
  All_Eye_Data <- rbind(All_Eye_Data, All_Eye_Data_VP)
  
}

write.csv2(All_Eye_Data,paste0(path, "Data/Tobii/All_Eye_Data.csv"))

# Percent invalid baseline 
count(All_Eye_Data[All_Eye_Data$blok==FALSE,])/count(All_Eye_Data[All_Eye_Data$blok,])

require(ez)

# Exclude bad baselines --> not exclude trials with bad baselines, but subjects with too many bad baselines 
table(All_Eye_Data$blok)
All_Eye_Data <- All_Eye_Data[All_Eye_Data$blok==TRUE,]
All_Eye_Data <- All_Eye_Data[All_Eye_Data$problem==0,]

##############################################
# Analyze aggregated data (across trials)
# 
#mean_All_Eye_Data <- aggregate(All_Eye_Data[,13:ncol(All_Eye_Data)], by=list(All_Eye_Data$vp, All_Eye_Data$cue), mean, na.rm=TRUE)
mean_All_Eye_Data <- All_Eye_Data %>%
  group_by(vp, cue)%>%
  select("1_CB":"8_FixDur")%>%
  summarise_all(mean, na.rm = TRUE)

write.csv2(mean_All_Eye_Data, paste0(path, "Data/Tobii/All_Eye_Data_Mean.csv"))  

##############################

# Eye-Tracking variables
cb <- mean_All_Eye_Data[,grep("CB",names(mean_All_Eye_Data))]
names(cb) <- paste("cb",1:ncol(cb),sep="")
cb <- data.frame(code=mean_All_Eye_Data[,1], cue=mean_All_Eye_Data[,2], cb)
cb <- pivot_longer(cb, cb1:cb8, names_to = "bin", values_to = "value")
cb$bin <- as.numeric(substr(cb$bin, 3, 3))
cb <- mutate (cb, 
              ActPass = case_when(
                cue == 1 | cue == 3 ~ "Passive",
                cue == 2 | cue == 4 ~ "Active"),
              Threat = case_when(
                cue == 1 | cue == 2 ~ "High",
                cue == 3 | cue == 4 ~ "Low"),
              cue = case_when(
                cue == 1 ~ "HighThreatPassive",
                cue == 2 ~ "HighThreatActive",
                cue == 3 ~ "LowThreatPassive",
                cue == 4 ~ "LowThreatActive"))
write.csv2(cb, paste0(path, "Data/Tobii/CB_Eye_Data_Mean.csv"), row.names=FALSE) 

fdur <- mean_All_Eye_Data[,grep("FixDur",names(mean_All_Eye_Data))]
names(fdur) <- paste("fixdur",1:ncol(fdur),sep="")
fdur <- data.frame(code=mean_All_Eye_Data[,1], cue=mean_All_Eye_Data[,2], fdur)
fdur <- pivot_longer(fdur, fixdur1:fixdur8, names_to = "bin", values_to = "value")
fdur$bin <- as.numeric(substr(fdur$bin, 7, 7))
fdur <- mutate (fdur, 
              ActPass = case_when(
                cue == 1 | cue == 3 ~ "Passive",
                cue == 2 | cue == 4 ~ "Active"),
              Threat = case_when(
                cue == 1 | cue == 2 ~ "High",
                cue == 3 | cue == 4 ~ "Low"),
              cue = case_when(
                cue == 1 ~ "HighThreatPassive",
                cue == 2 ~ "HighThreatActive",
                cue == 3 ~ "LowThreatPassive",
                cue == 4 ~ "LowThreatActive"))
write.csv2(fdur, paste0(path, "Data/Tobii/FixDur_Eye_Data_Mean.csv"), row.names=FALSE) 

fnum <- mean_All_Eye_Data[,grep("FixN",names(mean_All_Eye_Data))]
names(fnum) <- paste("fixnum",1:ncol(fnum),sep="")
fnum <- data.frame(code=mean_All_Eye_Data[,1], cue=mean_All_Eye_Data[,2], fnum)
fnum <- pivot_longer(fnum, fixnum1:fixnum8, names_to = "bin", values_to = "value")
fnum$bin <- as.numeric(substr(fnum$bin, 7, 7))
fnum <- mutate (fnum, 
              ActPass = case_when(
                cue == 1 | cue == 3 ~ "Passive",
                cue == 2 | cue == 4 ~ "Active"),
              Threat = case_when(
                cue == 1 | cue == 2 ~ "High",
                cue == 3 | cue == 4 ~ "Low"),
              cue = case_when(
                cue == 1 ~ "HighThreatPassive",
                cue == 2 ~ "HighThreatActive",
                cue == 3 ~ "LowThreatPassive",
                cue == 4 ~ "LowThreatActive"))
write.csv2(fnum, paste0(path, "Data/Tobii/FixNum_Eye_Data_Mean.csv"), row.names=FALSE) 

# Center Bias: 3 x 8 ANOVA (cue x second)------------------------------------------------------------------------

cb_anova = cb %>%
  mutate(bin = as.factor(bin), ActPass = as.factor(ActPass), Threat = as.factor(Threat))%>%
  ezANOVA(dv=.(value), wid=.(vp), within=.(bin,ActPass,Threat), detailed=T, type=3)

cb_anova %>% apa::anova_apa(force_sph_corr=T)

print(cb_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # bin
print(cb_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # ActPass
print(cb_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # Threat
print(cb_anova$ANOVA[5,] %>% partial_eta_squared_ci()) # bin * ActPass
print(cb_anova$ANOVA[6,] %>% partial_eta_squared_ci()) # bin * Threat
print(cb_anova$ANOVA[7,] %>% partial_eta_squared_ci()) # ActPass * Threat
print(cb_anova$ANOVA[8,] %>% partial_eta_squared_ci()) # bin * ActPass * Threat

cb_plot <- cb %>%
  group_by(bin, cue) %>%
  summarize(
    mean = mean(value),
    se = sd(value)/sqrt(n())
  )%>%
  ggplot(aes(x = bin, y = mean, color = cue, group = cue))+
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymax = mean+se, ymin = mean-se))
print(cb_plot)
ggsave(file = paste0(path,"cb_plot.png",sep=""), dpi = 600, width = 16.6, height = 16.6/2, units = "cm", bg="white",scale = 2)

# t-tests 
cb.mean = cb %>%
  group_by(vp, cue)%>%
  summarize(value = mean(value))
t.data.highshock= cb.mean %>%
  filter(cue == "HighThreatPassive")
t.data.flight = cb.mean %>%
  filter(cue == "HighThreatActive")
t.data.lowshock = cb.mean %>%
  filter(cue == "LowThreatPassive")
t.data.avoidance = cb.mean %>%
  filter(cue == "LowThreatActive")
#High shock vs low shock
t.test(t.data.highshock$value,t.data.lowshock$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)
#High shock vs. flight
t.test(t.data.highshock$value,t.data.flight$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)
#Low shock vs. avoidance
t.test(t.data.lowshock$value,t.data.avoidance$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)
#Flight vs. avoidance
t.test(t.data.avoidance$value,t.data.flight$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)

# Fixation Durations: 3 x 8 ANOVA (cue x second)-------------------------------------------------------------------

fdur_anova = fdur %>%
  mutate(bin = as.factor(bin), ActPass = as.factor(ActPass), Threat = as.factor(Threat))%>%
  ezANOVA(dv=.(value), wid=.(vp), within=.(bin,ActPass, Threat), detailed=T, type=3)

fdur_anova %>%
  apa::anova_apa(force_sph_corr=T)

print(fdur_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # bin
print(fdur_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # ActPass
print(fdur_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # Threat
print(fdur_anova$ANOVA[5,] %>% partial_eta_squared_ci()) # bin * ActPass
print(fdur_anova$ANOVA[6,] %>% partial_eta_squared_ci()) # bin * Threat
print(fdur_anova$ANOVA[7,] %>% partial_eta_squared_ci()) # ActPass * Threat
print(fdur_anova$ANOVA[8,] %>% partial_eta_squared_ci()) # bin * ActPass * Threat

fdur_plot <- fdur %>%
  group_by(bin, cue) %>%
  summarize(
    mean = mean(value),
    se = sd(value)/sqrt(n())
  )%>%
  ggplot(aes(x = bin, y = mean, color = cue, group = cue))+
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymax = mean+se, ymin = mean-se))
print(fdur_plot)

#means
fdur%>%
  group_by(ActPass, Threat)%>%
  summarize(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T))

fdur%>%
  group_by(vp, ActPass, Threat)%>%
  summarize(mean = mean(value, na.rm = T))%>%
  pivot_wider(names_from = Threat, values_from = mean)%>%
  mutate(diff = High - Low)%>%
  ungroup()%>%
  group_by(ActPass)%>%
  summarize(mean = mean(diff, na.rm = T),
            sd = sd(diff, na.rm = T))

# t-tests 
fdur.mean = fdur%>%
  group_by(vp, cue)%>%
  summarize(value = mean(value))
t.data.highshock= fdur.mean %>%
  filter(cue == "HighThreatPassive")
t.data.flight = fdur.mean %>%
  filter(cue == "HighThreatActive")
t.data.lowshock = fdur.mean %>%
  filter(cue == "LowThreatPassive")
t.data.avoidance = fdur.mean %>%
  filter(cue == "LowThreatActive")
#High shock vs low shock
t.test(t.data.highshock$value,t.data.lowshock$value,alternative = c("less"), paired=T) %>% apa::t_apa(es_ci=T)
#High shock vs. flight
t.test(t.data.highshock$value,t.data.flight$value,alternative = c("less"), paired=T) %>% apa::t_apa(es_ci=T)
#Low shock vs. avoidance
t.test(t.data.lowshock$value,t.data.avoidance$value,alternative = c("less"), paired=T) %>% apa::t_apa(es_ci=T)
#Flight vs. avoidance
t.test(t.data.avoidance$value,t.data.flight$value,alternative = c("less"), paired=T) %>% apa::t_apa(es_ci=T)

# Fixation Number: 3 x 8 ANOVA (cue x second) --------------------------------------------------------------------------

fnum_anova = fnum %>%
  mutate(bin = as.factor(bin), ActPass = as.factor(ActPass), Threat = as.factor(Threat))%>%
  ezANOVA(dv=.(value), wid=.(vp), within=.(bin,ActPass, Threat), detailed=T, type=3)

fnum_anova %>%
  apa::anova_apa(force_sph_corr=T)

print(fnum_anova$ANOVA[2,] %>% partial_eta_squared_ci()) # bin
print(fnum_anova$ANOVA[3,] %>% partial_eta_squared_ci()) # ActPass
print(fnum_anova$ANOVA[4,] %>% partial_eta_squared_ci()) # Threat
print(fnum_anova$ANOVA[5,] %>% partial_eta_squared_ci()) # bin * ActPass
print(fnum_anova$ANOVA[6,] %>% partial_eta_squared_ci()) # bin * Threat
print(fnum_anova$ANOVA[7,] %>% partial_eta_squared_ci()) # ActPass * Threat
print(fnum_anova$ANOVA[8,] %>% partial_eta_squared_ci()) # bin * ActPass * Threat

#means
fnum%>%
  group_by(vp, ActPass, Threat)%>%
  summarize(mean = mean(value, na.rm = T))%>%
  pivot_wider(names_from = Threat, values_from = mean)%>%
  mutate(diff = High - Low)%>%
  ungroup()%>%
  group_by(ActPass)%>%
  summarize(mean = mean(diff, na.rm = T),
            sd = sd(diff, na.rm = T))

fnum_plot <- fnum %>%
  group_by(bin, cue) %>%
  summarize(
    mean = mean(value),
    se = sd(value)/sqrt(n())
  )%>%
  ggplot(aes(x = bin, y = mean, color = cue, group = cue))+
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymax = mean+se, ymin = mean-se))
print(fnum_plot)

# t-tests 
fn.mean = fnum %>%
  group_by(vp, cue)%>%
  summarize(value = mean(value))
t.data.highshock= fn.mean %>%
  filter(cue == "HighThreatPassive")
t.data.flight = fn.mean %>%
  filter(cue == "HighThreatActive")
t.data.lowshock = fn.mean %>%
  filter(cue == "LowThreatPassive")
t.data.avoidance = fn.mean %>%
  filter(cue == "LowThreatActive")
#High shock vs low shock
t.test(t.data.highshock$value,t.data.lowshock$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)
#High shock vs. flight
t.test(t.data.highshock$value,t.data.flight$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)
#Low shock vs. avoidance
t.test(t.data.lowshock$value,t.data.avoidance$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)
#Flight vs. avoidance
t.test(t.data.avoidance$value,t.data.flight$value,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)

# # long to wide if necessary ------------------------------------------------------------------------------
# 
# cb_wide_out <- cb %>% 
#   group_by(vp,cue,bin) %>%  summarise(cb = mean(value), .groups = "drop") %>% 
#   mutate(bin = bin - 1.5, condition_names = paste0(cue,"_",bin)) %>% select(-cue, -bin) %>%
#   pivot_wider(., names_from = condition_names, values_from = cb)
# write.csv2(cb_wide_out,"CB_Daten.csv",row.names=F)
# 
# fixnum_wide_out <- fnum %>% 
#   group_by(vp,cue,bin) %>%  summarise(fnum = mean(value), .groups = "drop") %>% 
#   mutate(bin = bin - 1.5, condition_names = paste0(cue,"_",bin)) %>% select(-cue, -bin) %>%
#   pivot_wider(., names_from = condition_names, values_from = fnum)
# write.csv2(fixnum_wide_out,"FixNum_Daten.csv",row.names=F)
# 
# fixdur_wide_out <- fdur %>% 
#   group_by(vp,cue,bin) %>%  summarise(fdur = mean(value), .groups = "drop") %>% 
#   mutate(bin = bin - 1.5, condition_names = paste0(cue,"_",bin)) %>% select(-cue, -bin) %>%
#   pivot_wider(., names_from = condition_names, values_from = fdur)
# write.csv2(fixdur_wide_out,"FixDur_Daten.csv",row.names=F)


