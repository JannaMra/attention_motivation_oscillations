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

path <- paste0(getwd(),"/")
path.log <- "Data/log/" %>% paste0(path, .) 
path.seq <- "Data/sequences/" %>% paste0(path, .) 
savepath <- "Data/prot/" %>% paste0(path, .)  


vpn <- paste("vp",ifelse(1:52<10,"0",""),1:52,sep="")
# Exclusions:
exclusions <- c("vp05","vp22","vp30","vp51")
vpn <- vpn[!(vpn %in% c("vp05","vp22","vp30","vp51"))] #1 subject technical issues (underwent first block two times), 3 subjects did not follow the instructions (too many preliminary responses or too many responses in passive trials)


responses.summary <- data.frame()

for (vp in vpn) {
  #vp = "vp05"
  print(vp)
  
  seqdat1 <- read.table(paste(path.seq,vp,"_1.txt",sep=""),header=TRUE)
  seqdat2 <- read.table(paste(path.seq,vp,"_2.txt",sep=""),header=TRUE)
  seqdat <- bind_rows(seqdat1,seqdat2)
  seqdat$vp <- vp
  seqdat$trial <- c(1:160)
  
  seqdat$resp   <- NA  # Reaktion (Taste)
  seqdat$rt        <- NA  # Reaktionszeit
  seqdat$shock     <- 0   # Schock verabreicht?
  
  
  # Log einlesen
  header.logdat <- read.table(file=paste(path.log,vp,"_1-ShockScanning_Flight_New.log",sep=""), skip = 3, 
                              nrows = 1, header = FALSE, sep ="\t", stringsAsFactors = FALSE)
  if(vp == "vp26" | vp == "vp28"){
    logdat1  <- read.table(file=paste(path.log,vp,"_1-ShockScanning_Flight_New.log",sep=""), 
                           sep="\t", skip = 6, header=FALSE,fill=TRUE)
  } else{
    logdat1  <- read.table(file=paste(path.log,vp,"_1-ShockScanning_Flight_New.log",sep=""), 
                           sep="\t", skip = 4, header=FALSE,fill=TRUE) 
  }
  logdat2  <- read.table(file=paste(path.log,vp,"_2-ShockScanning_Flight_New.log",sep=""), 
                         sep="\t", skip = 9, header=FALSE,fill=TRUE)     # skip 8, so that "responses" are not counted as responses to cue
  logdat <- bind_rows(logdat1, logdat2) 
  header.logdat[3] <- "Event.Type"
  colnames( logdat ) <- unlist(header.logdat)
  logdat$Subject <- vp
  
  aktstim <- 1
  i <- 1
  
  while (i < nrow(logdat)) {
    if (as.character(logdat$Code[i])==paste(seqdat$pic[aktstim],".jpg",sep="")) {
      # Find next stimulus
      j <- i
      while ((j < nrow(logdat)) & (as.character(logdat$Code[j])!=paste(seqdat$pic[aktstim+1],".jpg",sep=""))) {
        j <- j + 1
      }

      # Find cue preceding stimulus
      cuepos <- i-1
      while ((cuepos > 1) & !(as.character(logdat$Code[cuepos]) %in% c("cue"))) {
        cuepos <- cuepos - 1
      }


      stimblock <- logdat[cuepos:(j-2),]

      # # Condition korrekt
      # if (stimblock$Code[1]!=c("flight", "shock", "noshock")[seqdat$cue[aktstim]+1]) {
      #   print("Wrong condition")
      # }

      # Shock given?
      seqdat$shock[aktstim] <- sum(stimblock$Code=="shock1")

      # Reaktionen finden und speichern
      flightprompt <- which(stimblock$Code=="cue")
      for (k in 1:nrow(stimblock)) {
        # Reaktion (nur die erste Reaktion werten)
        if ((stimblock$Event.Type[k]=="Response") & (is.na(seqdat$resp[aktstim]))) {
          reaktion <- as.numeric(as.character(stimblock$Code[k]))
          rtdat <-  stimblock$Time[k]-stimblock$Time[flightprompt]-100000 #stimblock$TTime[k]

          if (length(flightprompt)==1) {
            # Response occurred in flight trial
            seqdat$resp[aktstim]  <- reaktion
            seqdat$rt[aktstim] <- rtdat/10
          } else {
            # Response occurred in other trial
            seqdat$resp[aktstim]  <- 0
            rtdat <- stimblock$Time[k]-stimblock$Time[1]
            seqdat$rt[aktstim] <- rtdat/10
          }
        }
      }
    aktstim <- aktstim+1
    }
    i <- i+1
  }

  names(seqdat) <- c("pic","cue","iti","group", "vp", "trial", "resp","rt","shock")
  seqdat <- seqdat %>% select(vp, group, trial, pic, cue, iti, resp, rt, shock)
  seqdat$resp <- seqdat$resp %>% replace(is.na(.), 0)
  
  seqdat <- seqdat %>% 
       mutate(problem = case_when(
         cue==1 & resp==1 | cue==3 & resp == 1 | cue==2 & resp==0 | cue==4 & resp==0  | cue==2 & rt < 0 | cue == 4 & rt < 0 ~ 1,
         cue==1 & resp==0 | cue==3 & resp == 0 |cue==2 & resp==1 | cue==4 & resp==1 ~ 0)
       )
  write.csv2(seqdat,paste(savepath,vp,".csv",sep=""),row.names=FALSE,quote=FALSE)

  # Testberechnungen
  # Anzahl geschockte Trials per Condition
  print(table(seqdat$cue,seqdat$shock))
  # Anzahl Trials mit fehlerhaften Reaktionen
  print(paste0("No responses in flight trials: ",sum(seqdat$cue==2 & seqdat$resp==0)))
  print(paste0("No responses in avoidance trials: ",sum(seqdat$cue==4 & seqdat$resp==0)))
  print(paste0("Response in passive trials: ",sum((seqdat$cue==1 |seqdat$cue==3)& seqdat$resp==1)))
  # Anzahl Flight-Trials mit zu fruehen Reaktionen
  print(paste0("Premature flight responses: ",sum(seqdat$resp==1 & seqdat$rt<0,na.rm=TRUE)))

  No_flight_reponse <- sum(seqdat$cue==2 & seqdat$resp==0)
  Percent_successful_flight <- sum(seqdat$cue== 2 & seqdat$rt > 0 & seqdat$rt < 1000 & seqdat$shock==0)/40
  No_avoidance_response <- sum(seqdat$cue==4 & seqdat$resp==0)
  Percent_successful_avoidance <- sum(seqdat$cue== 4 & seqdat$rt > 0 & seqdat$rt < 1000 & seqdat$shock==0)/40
  Reaction_in_passive_Trial <- sum((seqdat$cue==1 | seqdat$cue==3) & seqdat$resp==1)
  Premature_responses <- sum((seqdat$cue==2 | seqdat$cue==4) & seqdat$resp==1 & seqdat$rt<0)
  RT_flight <- seqdat %>% 
    filter(cue == 2 & rt > 0 & rt < 1000) %>%
    summarize(
      mean = mean(rt))%>%
    as.double(mean)
  RT_avoidance <- seqdat %>% 
    filter(cue == 4 & rt > 0 & rt < 1000) %>%
    summarize(
      mean = mean(rt))%>%
    as.double(mean)
  summary.vp <- data.frame(vp, No_flight_reponse, Percent_successful_flight, No_avoidance_response, Percent_successful_avoidance, Reaction_in_passive_Trial, Premature_responses, RT_flight,  RT_avoidance )
  responses.summary <- bind_rows(responses.summary, summary.vp)
}

# mean 
responses.summary %>%
  summarise(flight_mean_rt = mean(RT_flight),
            sd = sd(RT_flight),
            flight_mean_success = mean(Percent_successful_flight),
            sd2 = sd(Percent_successful_flight),
            avoidance_mean_rt = mean(RT_avoidance),
            sd3 = sd(RT_avoidance),
            avoidance_mean_success = mean(Percent_successful_avoidance),
            sd4 = sd(Percent_successful_avoidance))

responses.summary.plot <- responses.summary %>%
  pivot_longer(cols = starts_with("RT"), names_to = "condition", names_prefix = "RT_", values_to = "rt")


response.plot <- ggplot(responses.summary.plot, aes(x=condition, y=rt, fill = condition)) + 
  geom_boxplot(color="black", notch=TRUE)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1,color="black", alpha = 0.3)+
  scale_fill_manual(values= alpha(c("darkblue", "#FF0000"), .6))+
  theme_classic()+
  labs(x = "Condition", y = "Mean Response Time")+
  scale_x_discrete(labels=c("Low Intensity", "High Intensity"))+
  theme(legend.position="none")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12))
response.plot

t.test(responses.summary$RT_flight,responses.summary$RT_avoidance,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)
t.test(responses.summary$Percent_successful_flight,responses.summary$Percent_successful_avoidance,alternative = c("greater"), paired=T) %>% apa::t_apa(es_ci=T)


for (vp in vpn) {
  #vp <- "vp01"
  cues <- read.csv2(paste(savepath,vp,".csv",sep=""))
  cues <- cues %>%
    mutate(trigger = ifelse(problem == 0, cue, 5))%>%
    select(trigger)
  write.table(cues,paste(savepath,vp,"triggerEEG.txt",sep=""),row.names=FALSE,quote=FALSE)
  
}
