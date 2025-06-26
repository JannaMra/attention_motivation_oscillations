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
library(afex)
library(apa)
library(lme4)
library(ggplot2)

path <- paste0(getwd(),"/")

trials.n = 160
sample.rate = 100
minhr =  45 #minimum plausible heart rate
maxhr = 120 #maximum plausible heart rate

scaling.window = c(-1, seq(0, 15, by=1)) # Scoring bins in seconds (real time scaling; may be non-integer)

scaleHR = function(hr_t, hr, st, en) {
  wsum = 0 #weighted sum (will be divided by valid ranges below)
  naRanges = 0
  i = which(hr_t > st)[1] #first r peak after marker start (might be NA)
  while (!is.na(i) && i <= length(hr_t) && hr_t[i] < en) { #build up weighted sum until last r peak in interval (exclusively! see below)
    #determine range
    lower = ifelse(i > 1 && hr_t[i-1] > st, hr_t[i-1], st) #find lower end of scoring interval (might be st if previous beat does not exist or out or scoring range)
    range = hr_t[i] - lower
    
    #determine value
    if (i > 1 && !is.na(hr[i])) value = hr[i]
    else {
      value = 0
      naRanges = naRanges + range #lower total range by invalid ranges
    }
    
    wsum = wsum + range * value
    i = i + 1
  } #wsum and naRanges are computed except for last beat
  
  #last beat
  if (!is.na(i) && i > 1 && i <= length(hr_t)) {
    range = en - max(c(hr_t[i-1], st)) #it might happen that there is only one beat before st and one after en => use hr[i] (implicitly done by setting range = en - st via max(c(...)))
    
    value = hr[i]
    if (is.na(value)) {
      value = 0
      naRanges = naRanges + range #lower total range by invalid ranges
    } else value = hr[i]
    
    wsum = wsum + range * value
  }
  
  result = wsum / (en - st - naRanges)
  if (is.na(result) || result ==  0) result = NA
  return(result)
}


# Read & Score HR --------------------------------------------------------------------
vpn.ecg.rpeaks = list.files(paste0(path,"Data/Physio/rpeaks/"), pattern=path.rpeaks.postfix, full.names=TRUE)
vpn.ecg.rpeaks <- vpn.ecg.rpeaks[!grepl("vp05|vp22|vp30|vp51", vpn.ecg.rpeaks)]  #filter general exclusions
#vpn.ecg.rpeaks <- vpn.ecg.rpeaks[!grepl("vp10|vp23", vpn.ecg.rpeaks)] #testweise Ausschluss falsche Anzahl an Trigger?? 
#ratings.all = read_rds("ratings.rds" %>% paste0(path.rds, .))
#rawfiles = vpn.ecg.rpeaks %>% gsub("rpeaks/", "", .) %>% gsub(path.rpeaks.postfix, "", .) %>% paste0(".txt")

hr.list = list() #vector("list", length(vpn.ecg.rpeaks))
for (vpi in seq(vpn.ecg.rpeaks)) {
  #vpi = 9 
  vp = vpn.ecg.rpeaks[vpi]
  #if (vp %>% pathToCode() %>% codeToNum() %in% exclusions.hr) next
  print(vp %>% pathToCode())
  
  code = vp %>% pathToCode() %>% pathToCode(file.ext = "_") %>% substr(1, 4)  
  #ratingfile = ratings.all %>% filter(subject == code %>% codeToNum())
  
  #load r peaks
  allrpeak = read.csv2(vp)[, 1] #only first column
  hr = 60/diff(allrpeak) #convert to bpm
  
  #load triggers
  trigger = read.table(paste0(path,"Data/Physio/EXP/",code,"_EKG.txt"),header = F) %>% mutate(trigger = ifelse(V3 %in% c(1,2,3,4),V3,0)) %>% .$trigger
  
  if (exclusions.phys.trials[[code]] %>% is.null() == F) { #exclude trials manually (cp. triggerCheck)
    triggers.only = trigger[trigger != 0]
    toExclude = exclusions.phys.trials[[code]]
    triggers.only[toExclude] = 0
    trigger[trigger != 0] = triggers.only
  }
  triggers.n = sum(trigger!=0)
  
  #exclude first wrong triggers for vp10 & vp23
  if(triggers.n > 160) {
    trigger[head(which(trigger > 0),triggers.n-160)] <- 0
    triggers.n = sum(trigger!=0)
  }
  
  conditions = trigger[trigger!=0]
  
  timeline = seq(trigger) / sample.rate
  marker = timeline[{trigger != 0} %>% which() %>% tail(trials.n)]
  
  missingBegin = min(marker) + min(scaling.window); missingBegin = missingBegin[missingBegin <= 0] # <= 0 because first sample point is 1 / sample.rate, i.e. 0 is out of range
  missingEnd = max(timeline) - (max(marker) + max(scaling.window)); missingEnd = missingEnd[missingEnd < 0] # < 0 because if difference exactly 0, then last sample point in rage
  
  #check if data is missing at the edges of data
  if (!is_empty(missingBegin)) warning(paste(code, ": Lacking data at BEGINNING", -missingBegin, "sec"))
  if (!is_empty(missingEnd)) warning(paste(code, ": Lacking data at END", -missingEnd, "sec"))
  
  #check for plausiblilty and issue warning
  if (min(hr) < minhr || max(hr) > maxhr) warning(paste(code, ": Implausible heart rate", hr %>% min() %>% round(1), "-", hr %>% max() %>% round(1)))
  
  hr = c(NA, hr) #matching allrpeak[i] to hr[i] (heart rate values only valid if PREVIOUS r peak exists)
  
  # Real time scoring
  allhr = numeric()
  for (trial in seq(marker)) {
    #trial = 1
    mtime = marker[trial]
    hr_t = allrpeak - mtime #heart rate time (relative to marker)
    
    # Real time scaling
    hrtrial = numeric()
    for (j in seq(scaling.window)[-1]) { #for all indices except for the first (due to j-1 indexing)
      current = ifelse(mtime + scaling.window[j] < 0, NA, #skip marker time points that refer to negative times (i.e. out of data)
                       scaleHR(hr_t, hr, scaling.window[j-1], scaling.window[j]))
      hrtrial = c(hrtrial, current)
    }
    
    allhr = rbind(allhr, hrtrial)
  }
  
  deltaval = allhr - matrix(allhr[, 1], nrow=nrow(allhr), ncol=ncol(allhr))
  #ratings.conditions = ratingfile$condition %>% tail(marker %>% length())
  #shocks = ratingfile$shock == "True" %>% tail(n=trials.n)
  out = data.frame(trial = 1:trials.n, condition = conditions, 
                   #shock = shocks, shockPrior = c(FALSE, lag(shocks)[-1]),
                   hrbl = allhr[, 1], hr = deltaval[, 2:ncol(deltaval)])
  
  hr.list[[code]] = out
  #write.csv2(out, paste(savepath.ecg, code,"_task.csv",sep=""), row.names=FALSE, quote=FALSE)
}

#list to one giant dataframe
heart.wide = hr.list %>% bind_rows(.id="subject") %>% 
  mutate(subject = subject %>% gsub("\\D", "", .) %>% as.integer()) %>% tibble() %>% 
  group_by(subject, condition) %>% mutate(trial_condition = 1:n()) %>% ungroup() %>% select(subject, trial, condition, trial_condition, everything())
rm(hr.list); row.names(heart.wide) = NULL

heart = heart.wide %>% gather(key="time", value="HRchange", matches("hr\\.\\d+")) %>% tibble() %>% 
  mutate(time = time %>% gsub("hr.", "", .) %>% as.integer() %>% as.factor(),
         condition = as.factor(condition),
         ActPass = case_when(condition ==1 | condition == 3 ~ "Passive",
                           condition == 2 | condition == 4 ~ "Active"),
       Intensity = case_when(condition == 1 | condition == 2 ~ "High",
                             condition== 3 | condition == 4 ~ "Low"))%>%
  mutate(ID = paste0("vp",str_pad(subject,2,pad="0")))

heart.data <- heart.wide %>% mutate(ID = paste0("vp",str_pad(subject,2,pad="0"))) %>% select(ID, trial, hrbl:hr.15)

#write_rds(heart, "Heart_df.rds")

#######################################################
## FOR DATA ANALYSIS READ PREPROCESSED FILES HERE #####
#######################################################

heart = read_rds("heart_df.rds" )

# Inference statistics--------------------------------------------------------------------

# read in prot files 
All_Prot_Data <- read.csv2(paste0(path.prot,"All_Prots.csv"))

# filter out problematic trials 

heart <- heart %>% rename(vp = ID)
heart <- full_join(All_Prot_Data, heart)

heart <- heart %>%
  filter(problem == 0)


anova <- heart %>%  
  filter(!(time %in% c(11,12,13,14,15))) %>%
  #filter(!(time %in% c(10.5,11,11.5,12,12.5,13,13.5,14,14.5,15))) %>%
  group_by(subject,time,ActPass, Intensity) %>% 
  summarize(HRchange = mean(HRchange, na.rm = T)) %>% ungroup() %>%
  mutate(subject = as.factor(subject), ActPass = as.factor(ActPass), Intensity = as.factor(Intensity), time = as.factor(time)) %>%
  ez::ezANOVA(dv=.(HRchange), wid=.(subject), 
              within=.(ActPass, Intensity, time), 
              #between=.(pairs),
              detailed=T, type=2)


anova %>% anova_apa()
anova$ANOVA[2,] %>% partial_eta_squared_ci() #ActPass
anova$ANOVA[3,] %>% partial_eta_squared_ci() #Intensity
anova$ANOVA[4,] %>% partial_eta_squared_ci() #time
anova$ANOVA[5,] %>% partial_eta_squared_ci() #ActPass * Intensity
anova$ANOVA[6,] %>% partial_eta_squared_ci() #ActPass * time
anova$ANOVA[7,] %>% partial_eta_squared_ci() #Intensity * time
anova$ANOVA[8,] %>% partial_eta_squared_ci() #ActPass * Intensity * time

# t-tests 
heart.mean = heart %>%  
  filter(!(time %in% c(1,2,3,4,5,6,11,12,13,14,15))) %>%
  #filter(!(time %in% c(10.5,11,11.5,12,12.5,13,13.5,14,14.5,15))) %>%
  group_by(subject,ActPass, Intensity) %>% 
  summarize(value = mean(HRchange, na.rm = T))
t.data.passive.high= heart.mean %>%
  filter(ActPass == "Passive" & Intensity == "High")
t.data.passive.low= heart.mean %>%
  filter(ActPass == "Passive" & Intensity == "Low")
#High shock vs low shock
t.test(t.data.passive.high$value,t.data.passive.low$value,alternative = c("less"), paired=T) %>% apa::t_apa(es_ci=T)


hr_wide_out <- heart %>%
  #filter(ID %in% responder) %>%
  select(ID, condition, time, HRchange) %>%
  group_by(ID,condition,time) %>%  summarise(HRchange = mean(HRchange, na.rm = T), .groups = "drop") %>%
  mutate(condition = factor(condition, levels=c(1,2,3,4),labels =c("HI_threat","HI_flight","LI_threat", "LI_flight")),
         condition_names = paste0(condition,"_",time)) %>% select(-condition, -time) %>%
  pivot_wider(., names_from = condition_names, values_from = HRchange)

write.csv2(hr_wide_out,"HR_Daten.csv",row.names=F)


# heart %>%  
#   filter(!(time %in% c(10.5,11,11.5,12,12.5,13,13.5,14,14.5,15))) %>%
#   group_by(subject,time,condition) %>% 
#   summarize(HRchange = mean(HRchange)) %>%
#   mutate(subject = as.factor(subject), condition = as.factor(condition), time = as.factor(time)) %>%
#   aov_ez(dv="HRchange", id="subject", . ,within = c("condition","time")) %>% 
#   emmeans::emmeans(~ condition|time) %>%
#   pairs() %>% summary(adjust = "FDR")





# heart %>% 
#   filter(!(time %in% c(11,11.5,12,12.5,13,13.5,14,14.5,15))) %>%
#   group_by(condition, time) %>% 
#   summarise(HRchange.se = sd(HRchange, na.rm=T)/sqrt(n()), HRchange = mean(HRchange, na.rm=T)) %>% 
#   mutate(time = time %>% as.character() %>% as.numeric()) %>% 
#   #bind_rows(data.frame(condition=unique(heart$condition), time=0, HRchange=0, HRchange.se=0)) %>% #add origin
#   ggplot(aes(x=time-0.5,y=HRchange, group = condition, color = condition)) +
#   #geom_rect(aes(xmin=4, xmax=6, ymin=-2, ymax=1), fill = "deepskyblue1", alpha=0.03) +
#   #geom_errorbar(aes(ymin=HRchange-HRchange.se,ymax=HRchange+HRchange.se), width=0.2)+
#   geom_path()+
#   #geom_point(size=3, shape=21)+
#   scale_color_manual(values=c("red","orange","darkgreen"), label = c("Threat","Flight","Safety")) +
#   scale_x_continuous(name ="Time [s]", limits = c(0,10), breaks = c(0,1,2,3,4,5,6,7,8,9,10)) + 
#   scale_y_continuous(name="Heart rate change [bpm]")+
#   theme_classic() +
#   theme(legend.position=c(0.15,0.18),
#         legend.title = element_blank(),
#         #panel.grid.major.y = element_line(size=0.5, color="#DDDDDD")
#   )

# ggsave(paste0("../Plots/HR/acq_cs_hr_bl_lines.png"), type="cairo-png", width=1920/400, height=1080/300, dpi=300)
# 
# # with US 
# heart %>% 
#   #filter(!(time %in% c(11,11.5,12,12.5,13,13.5,14,14.5,15))) %>%
#   group_by(condition, time) %>% 
#   summarise(HRchange.se = sd(HRchange, na.rm=T)/sqrt(n()), HRchange = mean(HRchange, na.rm=T)) %>% 
#   mutate(time = time %>% as.character() %>% as.numeric()) %>% 
#   #bind_rows(data.frame(condition=unique(heart$condition), time=0, HRchange=0, HRchange.se=0)) %>% #add origin
#   ggplot(aes(x=time-0.5,y=HRchange, group = condition, color = condition)) +
#   #geom_rect(aes(xmin=4, xmax=6, ymin=-2, ymax=1), fill = "deepskyblue1", alpha=0.03) +
#   geom_ribbon(aes(ymin=HRchange-HRchange.se,ymax=HRchange+HRchange.se, fill = condition), color = "transparent", alpha=0.2, show_guide = FALSE)+
#   geom_path()+
#   #geom_point(size=3, shape=21)+
#   scale_color_manual(values=c("red","orange","darkgreen"), label = c("Threat","Flight","Safety")) +
#   scale_fill_manual(values=c("red","orange","darkgreen"), label = c("Threat","Flight","Safety")) +
#   scale_x_continuous(name ="Time [s]") + 
#   scale_y_continuous(name="Heart rate change [bpm]")+
#   theme_classic() +
#   theme(legend.position=c(0.15,0.18),
#         legend.title = element_blank(),
#         #panel.grid.major.y = element_line(size=0.5, color="#DDDDDD")
#   )
# 
# ggsave(paste0("../Plots/HR/acq_cs_hr_bl_lines.png"), type="cairo-png", width=1920/400, height=1080/300, dpi=300)
# 


# 
# 
# 
# 
# hr_out <- heart %>% group_by(subject, condition, time) %>%
#   summarise(hr = mean(HRchange)) %>%
#   mutate(time = as.numeric(time),
#          ID = code, #paste0("tfo",str_pad(subject,2,pad="0")),
#          condition = factor(condition, levels = c(1,2,3), labels = c("threat", "flight", "safety")),
#          names = paste0(condition,"_",time)) %>% select(-subject) %>%
#   pivot_wider(id_cols = ID, names_from = names, values_from = hr) 
# hr_out <- hr_out[order(hr_out$ID),]
# write.csv2(hr_out,"HR_Daten.csv", row.names = F)
# 
# 
# 
# heart_wue <- heart %>% group_by(subject, condition, time) %>%
#   summarise(hr = mean(HRchange)) %>% ungroup() %>%
#   mutate(bin = as.numeric(time),
#          ID = subject,
#          cue = factor(condition, levels = c(1,2,3), labels = c("threat", "flight", "safety"))) %>% 
#   select(ID, cue, bin, hr, -subject) %>% mutate(site = "wuerzburg")
# 
# 
# # Gainesville Data
# 
# library(R.matlab)
# 
# hr_gainesville <- readMat("../Florida Sample/HR/allsubjects.mat")
# 
# hr_gainesville1 <- data.frame(hr_gainesville[[1]]) %>% 
#   mutate(code = c(1:49,51,52),
#          cue = "threat") %>%
#   rowwise() %>%
#   mutate(t0 = X1, t1 = mean(X2,X3), t2 = mean(X4,X5),t3 = mean(X6,X7),t4 = mean(X8,X9),
#          t5 = mean(X10,X11),t6 = mean(X12,X13),t7 = mean(X14,X15),t8 = mean(X16,X17), 
#          t9 = mean(X18,X19),t10 = mean(X20,X21)) %>%
#   select(-(X1:X22)) %>%
#   pivot_longer(t0:t10, names_to = "bin", values_to = "hr") %>%
#   mutate(bin = as.numeric(str_remove(bin,"t")))
# 
# hr_gainesville2 <- data.frame(hr_gainesville[[2]]) %>% 
#   mutate(code = c(1:49,51,52),
#          cue = "flight") %>%
#   pivot_longer(X1:X22, names_to = "bin", values_to = "hr") %>%
#   mutate(bin = as.numeric(str_remove(bin,"X"))-1)
# 
# 
# hr_gainesville2 <- data.frame(hr_gainesville[[2]]) %>% 
#   mutate(code = c(1:49,51,52),
#          cue = "flight") %>%
#   rowwise() %>%
#   mutate(t0 = X1, t1 = mean(X2,X3), t2 = mean(X4,X5),t3 = mean(X6,X7),t4 = mean(X8,X9),
#          t5 = mean(X10,X11),t6 = mean(X12,X13),t7 = mean(X14,X15),t8 = mean(X16,X17), 
#          t9 = mean(X18,X19),t10 = mean(X20,X21)) %>%
#   select(-(X1:X22)) %>%
#   pivot_longer(t0:t10, names_to = "bin", values_to = "hr") %>%
#   mutate(bin = as.numeric(str_remove(bin,"t")))
# 
# 
# hr_gainesville3 <- data.frame(hr_gainesville[[3]]) %>% 
#   mutate(code = c(1:49,51,52),
#          cue = "safety") %>%
#   pivot_longer(X1:X22, names_to = "bin", values_to = "hr") %>%
#   mutate(bin = as.numeric(str_remove(bin,"X"))-1)
# 
# 
# hr_gainesville3 <- data.frame(hr_gainesville[[3]]) %>% 
#   mutate(code = c(1:49,51,52),
#          cue = "safety") %>%
#   rowwise() %>%
#   mutate(t0 = X1, t1 = mean(X2,X3), t2 = mean(X4,X5),t3 = mean(X6,X7),t4 = mean(X8,X9),
#          t5 = mean(X10,X11),t6 = mean(X12,X13),t7 = mean(X14,X15),t8 = mean(X16,X17), 
#          t9 = mean(X18,X19),t10 = mean(X20,X21)) %>%
#   select(-(X1:X22)) %>%
#   pivot_longer(t0:t10, names_to = "bin", values_to = "hr") %>%
#   mutate(bin = as.numeric(str_remove(bin,"t")))
# 
# 
# hr_ga <- rbind(hr_gainesville1,hr_gainesville2) %>% rbind(.,hr_gainesville3) %>%
#   mutate(code = code+100) %>% rename(ID = code)
# 
# 
# 
# 
# 
# hr_ga %>% 
#   mutate(cue = factor(cue, levels = c("threat","flight","safety"), labels = c("threat", "flight", "safety"))) %>%
#   group_by(cue, bin) %>% 
#   summarise(hr.se = sd(hr, na.rm=T)/sqrt(n()), hr = mean(hr, na.rm=T)) %>% 
#   #bind_rows(data.frame(cue=unique(heart$cue), time=0, HRchange=0, HRchange.se=0)) %>% #add origin
#   ggplot(aes(x=bin,y=hr, group = cue, color = cue)) +
#   #geom_rect(aes(xmin=4, xmax=6, ymin=-2, ymax=1), fill = "deepskyblue1", alpha=0.03) +
#   geom_ribbon(aes(ymin=hr-hr.se,ymax=hr+hr.se, fill = cue), color = "transparent", alpha=0.1, show_guide = FALSE)+
#   geom_path()+
#   #geom_point(size=3, shape=21)+
#   scale_color_manual(values=c("red","orange","darkgreen"), label = c("Threat","Flight","Safety")) +
#   scale_fill_manual(values=c("red","orange","darkgreen"), label = c("Threat","Flight","Safety")) +
#   #scale_x_continuous(name ="Time [s]", breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)) + 
#   scale_y_continuous(name="Heart rate change [bpm]")+
#   theme_classic() +
#   theme(legend.position=c(0.15,0.88),
#         legend.title = element_blank(),
#         #panel.grid.major.y = element_line(size=0.5, color="#DDDDDD")
#   )
# 
# 
# anova <- hr_ga %>%  
#   filter(!(bin %in% c(0))) %>%
#  # filter((bin %in% c(9,10))) %>%
#   group_by(ID,bin,cue) %>% 
#   summarize(hr = mean(hr)) %>% ungroup() %>%
#   mutate(subject = as.factor(subject), cue = as.factor(cue), bin = as.factor(bin)) %>%
#   ez::ezANOVA(dv=.(hr), wid=.(ID), 
#               within=.(cue, bin), 
#               #between=.(pairs),
#               detailed=T, type=2)
# 
# 
# anova %>% anova_apa()
# anova$ANOVA[2,] %>% partial_eta_squared_ci()
# anova$ANOVA[3,] %>% partial_eta_squared_ci()
# anova$ANOVA[4,] %>% partial_eta_squared_ci()
# 
# 
# # All together
# 
# hr_all <- rbind(heart_wue, hr_ga %>% filter(time!=0) %>% mutate(site = "gainesville"))
# 
# 
# 
# hr_all %>% 
#   filter(time < 21) %>% 
#   group_by(cue, bin) %>% 
#   summarise(hr.se = sd(hr, na.rm=T)/sqrt(n()), hr = mean(hr, na.rm=T)) %>% 
#   #bind_rows(data.frame(condition=unique(heart$condition), time=0, HRchange=0, HRchange.se=0)) %>% #add origin
#   ggplot(aes(x=bin/2,y=hr, group = cue, color = cue)) +
#   #geom_rect(aes(xmin=4, xmax=6, ymin=-2, ymax=1), fill = "deepskyblue1", alpha=0.03) +
#   geom_errorbar(aes(ymin=hr-hr.se,ymax=hr+hr.se), width=0.2)+
#   geom_path()+
#   #geom_point(size=3, shape=21)+
#   scale_color_manual(values=c("darkred","orange","blue")) +
#   #scale_x_continuous(name ="Time [s]", breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)) + 
#   scale_y_continuous(name="Heart rate change [bpm]")+
#   theme_classic() +
#   theme(legend.position=c(0.15,0.18),
#         legend.title = element_blank(),
#         #panel.grid.major.y = element_line(size=0.5, color="#DDDDDD")
#   )
# 
# 
# anova <- hr_all %>%  
#   filter(bin != 0) %>%
#   filter(bin < 11) %>%
#   group_by(ID,bin,cue) %>% 
#   summarize(hr = mean(hr)) %>% ungroup() %>%
#   mutate(ID = as.factor(ID), cue = as.factor(cue), bin = as.factor(bin)) %>%
#   ez::ezANOVA(dv=.(hr), wid=.(ID), 
#               within=.(cue, bin), 
#               #between=.(pairs),
#               detailed=T, type=2)
# 
# 
# anova %>% anova_apa()
# anova$ANOVA[2,] %>% partial_eta_squared_ci()
# anova$ANOVA[3,] %>% partial_eta_squared_ci()
# anova$ANOVA[4,] %>% partial_eta_squared_ci()
