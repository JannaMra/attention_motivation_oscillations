
library(tidyverse)
library(cowplot)
library(ggpubr)
library(svglite)
library(patchwork)
library(ggplot2)

# paths
path <- paste0(getwd(),"/")

# Add lines indicating pairwise significant differences

add_cluster_analysis.actpass <- function(lineplot, clusters_within, comparison, col_ph="black", binwidth=0) {
  for (i in 1:nrow(clusters_within)) {
    if (clusters_within$p[i] < 0.05 & clusters_within$cue[i]==comparison) {
      x <- clusters_within$start[i]-binwidth; xend <- clusters_within$end[i]
      if (grepl("high", comparison)) {
        y <- y1 ; yend <- y1} else {
          y <- y4 ; yend <- y4
        }
      lineplot <- lineplot +
        geom_segment(mapping=aes_string(x=x, y=y, xend=xend, yend=yend),colour=col_ph,size=1)
    }
  }
  return(lineplot)
}

add_cluster_analysis.highlow <- function(lineplot, clusters_within, comparison, col_ph="black", binwidth=0) {
  for (i in 1:nrow(clusters_within)) {
    if (clusters_within$p[i] < 0.05 & clusters_within$cue[i]==comparison) {
      x <- clusters_within$start[i]-binwidth; xend <- clusters_within$end[i]
      if (grepl("active", comparison)) {
        y <- y2 ; yend <- y2
        lineplot <- lineplot +
          geom_segment(mapping=aes_string(x=x, y=y, xend=xend, yend=yend),colour=col_ph, linetype = "dotdash",size=1)} else {
          y <- y3 ; yend <- y3
          lineplot <- lineplot +
            geom_segment(mapping=aes_string(x=x, y=y, xend=xend, yend=yend),colour=col_ph,size=1)
        }
    }
  }
  return(lineplot)
}


## Prepare autonomic data

eda_long <-  readRDS("EDA_long.RData") #read.table("EDA_Daten.csv", sep=";", dec=",", header = TRUE, check.names=F)
sc_long <- eda_long
names(sc_long) <- c("vp","cue","bin","value", "ActPass", "Intensity")
sc_long <- sc_long %>%
 group_by(vp, cue, bin) %>% summarise(value = mean(value))
sc_long$bin <- sc_long$bin - 0.5 
#sc <- data.frame(pivot_wider(sc_long, names_from="bin", values_from = "value"))

# Design:
# cue 1: high intensity shock (shock)
# cue 2: high intensity flight (flight)
# cue 3: low intensity shock (safety)
# cue 4: low intensity flight (active control)

sc_plot <- 
  sc_long %>%
  group_by(bin,cue) %>%
  summarise(mean = mean(value, na.rm= TRUE), se = sd(value,na.rm=TRUE) / sqrt(n())) %>%
  ggplot(aes(x = bin, group = cue, y = mean)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE")+
  geom_line(aes(colour = cue, linetype = cue), linewidth = 0.8) +
  scale_linetype_manual(name="Trial type", values=c("solid", "dotdash", "solid", "dotdash"), labels = c("HI Shock","HI Flight","LI Shock", "LI Flight"))+
  scale_color_manual(name="Trial type", values = c("#FF0000","#FF0000","#000099", "#000099"),labels = c("HI Shock","HI Flight","LI Shock", "LI Flight")) + 
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = cue), alpha = 0.1) + 
  scale_fill_manual(values = c("#FF0000","#FF0000","#000099", "#000099")) +
  guides(fill="none") +
  ylab(bquote('Δ SC (μS)')) +
  scale_x_continuous("Time (s)",limits= c(0,15),breaks=c(0,2,4,6,8,10,12,14)) +
  scale_y_continuous(limits= c(-0.45,0.2)) +
  ggtitle("Skin conductance")+
  theme_classic()+
  theme(legend.position = "none")

y1 <- ggplot_build(sc_plot)$layout$panel_params[[1]]$y$limits[2]
y2 <- y1 - diff(ggplot_build(sc_plot)$layout$panel_params[[1]]$y$limits)/25
y3 <- y2 - diff(ggplot_build(sc_plot)$layout$panel_params[[1]]$y$limits)/25
y4 <- y3 - diff(ggplot_build(sc_plot)$layout$panel_params[[1]]$y$limits)/25

# Add results of post-hoc tests
sc_plot <- add_cluster_analysis.actpass(sc_plot, eda_clusters_within, "eda.high.act_vs_pass", "#FF0000", 1)
sc_plot <- add_cluster_analysis.actpass(sc_plot, eda_clusters_within, "eda.low.act_vs_pass", "#000099", 1)
sc_plot <- add_cluster_analysis.highlow(sc_plot, eda_clusters_within, "eda.passive.high_vs_low", "#b74ac0", 1)
sc_plot <- add_cluster_analysis.highlow(sc_plot, eda_clusters_within, "eda.active.high_vs_low", "#b74ac0", 1)
# sc_plot <- sc_plot +
#   annotate("text", x=0, y=y1, label= "HI active vs. HI passive", size = 4, color = "#FF0000", hjust = 0)+
#   annotate("text", x=0, y=y2, label= "HI active vs. LI active", size = 4, color = "#b74ac0", hjust = 0)+
#   annotate("text", x=0, y=y3, label= "HI passive vs. LI passive", size = 4, color = "#b74ac0", hjust = 0)+
#   annotate("text", x=0, y=y4, label= "LI active vs. LI passive", size = 4, color = "#000099", hjust = 0)

#ggsave(file = paste(savepath,"SC_timecourse_faces.svg",sep=""), dpi = 600, width = 10, height = 10, units = "in", bg="white")
#ggsave(file = paste(savepath,"SC_timecourse_faces.pdf",sep=""), dpi = 600, width = 10, height = 10, units = "in", bg="white")

sc_plot 


#same for HR data
#hr <-  read.table("HR_Daten.csv", sep=";", dec=",", header = TRUE, check.names=F)
hr <- read_rds("heart_df.rds")%>%
  rename(code = vp, bin = time)%>%
  mutate(cue=as.factor(cue))%>%
  group_by(code,bin,cue, ActPass, Intensity)%>%
  summarize(
    value = mean(HRchange, na.rm = T)
  )%>% ungroup()
#hr_long <- pivot_longer(hr, cols = threat_1:safety_30, names_sep = "_", names_to = c("cue", "bin"), values_to = "value" ) %>%
  #mutate(bin = as.numeric(bin)) %>% filter(bin <= 20)



#names(hr) <- c("code","cue","bin","value")

hr$bin <- as.numeric(hr$bin) - 0.5
#hr <- data.frame(pivot_wider(hr_long, names_from="bin", values_from = "value"))

hr_plot <- hr %>%
  group_by(bin,cue) %>%
  summarise(mean = mean(value, na.rm= TRUE), se = sd(value,na.rm=TRUE) / sqrt(n())) %>%
  ggplot(aes(x = bin, group = cue, y = mean)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE")+
  geom_line(aes(colour = cue, linetype = cue), linewidth = 0.8) +
  scale_linetype_manual(name="Trial type", values=c("solid", "dotdash", "solid", "dotdash"), labels = c("HI Shock","HI Flight","LI Shock", "LI Flight"))+
  scale_color_manual(name="Trial type", values = c("#FF0000","#FF0000","#000099", "#000099"),labels = c("HI Shock","HI Flight","LI Shock", "LI Flight")) + 
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = cue), alpha = 0.2) +
  scale_fill_manual(values = c("#FF0000","#FF0000","#000099", "#000099")) +
  guides(fill="none") +
  ylab(bquote('Δ HR (bpm)')) +
  xlab("Time (s)") +
  scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14)) +
  scale_y_continuous(limits = c(-5,4)) +
  ggtitle("Heart rate") +
  theme_classic()+
  theme(legend.position = "none")

#y1 <- ggplot_build(hr_plot)$layout$panel_params[[1]]$y$limits[2]
#y2 <- y1 - diff(ggplot_build(hr_plot)$layout$panel_params[[1]]$y$limits)/50

# hr_p <- rbind(condcomp(hr[hr$cue=="flight",3:ncol(hr)],hr[hr$cue=="threat",3:ncol(hr)],"fdr"),
#               condcomp(hr[hr$cue=="flight",3:ncol(hr)],hr[hr$cue=="safety",3:ncol(hr)],"fdr"))
# hr_ph <- data.frame(cond=rep(c("fsh","fsa"),each=ncol(hr_p)),
#                     x=seq(0.5,15,1),y=rep(c(y1,y2),each=ncol(hr_p)),
#                     psig=c(hr_p[1,]<0.05,hr_p[2,]<0.05))

# # Add results of post-hoc tests
# hr_plot <- add_phtests(hr_plot, hr_ph, "fsh", "red1", 0.25)
# hr_plot <- add_phtests(hr_plot, hr_ph, "fsa", "limegreen", 0.25)

y1 <- ggplot_build(hr_plot)$layout$panel_params[[1]]$y$limits[2]
y2 <- y1 - diff(ggplot_build(hr_plot)$layout$panel_params[[1]]$y$limits)/25
y3 <- y2 - diff(ggplot_build(hr_plot)$layout$panel_params[[1]]$y$limits)/25
y4 <- y3 - diff(ggplot_build(hr_plot)$layout$panel_params[[1]]$y$limits)/25

# Add results of post-hoc tests
hr_plot <- add_cluster_analysis.actpass(hr_plot, hr_clusters_within, "hr.high.act_vs_pass", "#FF0000", 1)
hr_plot <- add_cluster_analysis.actpass(hr_plot, hr_clusters_within, "hr.low.act_vs_pass", "#000099", 1)
hr_plot <- add_cluster_analysis.highlow(hr_plot, hr_clusters_within, "hr.passive.high_vs_low", "#b74ac0", 1)
hr_plot <- add_cluster_analysis.highlow(hr_plot, hr_clusters_within, "hr.active.high_vs_low", "#b74ac0", 1)

hr_plot



#Pupil data
pupil_long <-  read.csv2(paste0(path.pupil, "pupil_long.csv"))%>% select(-1)%>%
  mutate(cue = as.factor(cue))

pupil_long$bin <- as.numeric(pupil_long$bin)-0.25
#pupil_long$bin <- (pupil_long$bin)/2 - 0.25 

pupil_plot <- pupil_long %>%
  group_by(bin,cue) %>%
  summarise(mean = mean(value, na.rm= TRUE), se = sd(value,na.rm=TRUE) / sqrt(n())) %>%
  ggplot(aes(x = bin, group = cue, y = mean)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE")+
  geom_line(aes(colour = cue, linetype = cue), linewidth = 0.8) +
  scale_linetype_manual(name="Trial type", values=c("solid", "dotdash", "solid", "dotdash"), labels = c("HI Shock","HI Flight","LI Shock", "LI Flight"))+
  scale_color_manual(name="Trial type", values = c("#FF0000","#FF0000","#000099", "#000099"),labels = c("HI Shock","HI Flight","LI Shock", "LI Flight")) + 
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = cue), alpha = 0.1) + 
  scale_fill_manual(values = c("#FF0000","#FF0000","#000099", "#000099")) +
  guides(fill="none") +
  ylab(bquote('Δ Pupil diameter (mm)')) +
  xlab("Time (s)") +  
  scale_x_continuous("Time (s)",limits= c(0,10),breaks=c(0,2,4,6,8,10)) +
  scale_y_continuous(limits=c(-0.4,0.15)) +
  ggtitle("Pupil dilation")+
  theme_classic()+
  theme(legend.position = "none")

y1 <- ggplot_build(pupil_plot)$layout$panel_params[[1]]$y$limits[2]
y2 <- y1 - diff(ggplot_build(pupil_plot)$layout$panel_params[[1]]$y$limits)/25
y3 <- y2 - diff(ggplot_build(pupil_plot)$layout$panel_params[[1]]$y$limits)/25
y4 <- y3 - diff(ggplot_build(pupil_plot)$layout$panel_params[[1]]$y$limits)/25

# Add results of post-hoc tests
pupil_plot <- add_cluster_analysis.actpass(pupil_plot, pupil_clusters_within, "pupil.high.act_vs_pass", "#FF0000", 0.5)
pupil_plot <- add_cluster_analysis.actpass(pupil_plot, pupil_clusters_within, "pupil.low.act_vs_pass", "#000099", 0.5)
pupil_plot <- add_cluster_analysis.highlow(pupil_plot, pupil_clusters_within, "pupil.passive.high_vs_low", "#b74ac0", 0.5)
pupil_plot <- add_cluster_analysis.highlow(pupil_plot, pupil_clusters_within, "pupil.active.high_vs_low", "#b74ac0", 0.5)


pupil_plot

#all_physio <- hr_plot + sc_plot/pupil_plot + plot_layout(widths = c(2,1))
#all_physio

# Prepare Fixation Data & Center Bias

cb <- read.csv2(paste0(path, "Data/Tobii/CB_Eye_Data_Mean.csv"))
fdur <- read.csv2(paste0(path, "Data/Tobii/FixDur_Eye_Data_Mean.csv"))
fnum <- read.csv2(paste0(path, "Data/Tobii/FixNum_Eye_Data_Mean.csv"))

# increase bin by 2s for each value for correct visual representation in plot
fdur$bin <- fdur$bin + 1.5 
fnum$bin <- fnum$bin + 1.5
cb$bin <- cb$bin + 1.5

# Plot Number of Fixations

fixnum_plot <- fnum %>%
  group_by(bin,cue) %>%
  summarise(
    mean = mean(value, na.rm= TRUE), 
    se = sd(value,na.rm=TRUE) / sqrt(n())) %>%
  ggplot(aes(x = bin, group = cue, y = mean)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE")+
  geom_line(aes(colour = cue, linetype = cue), size = 0.8) +
  scale_linetype_manual(name="Trial type", values=c( "dotdash","solid",  "dotdash", "solid"), labels = c("HI active","HI passive", "LI active","LI passive"))+
  scale_color_manual(name="Trial type", values = c("#FF0000","#FF0000","#000099", "#000099"),labels = c("HI active","HI passive", "LI active","LI passive")) + 
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = cue), alpha = 0.1) + 
  scale_fill_manual(values = c("#FF0000","#FF0000","#000099", "#000099")) +
  guides(fill=FALSE) +
  ylab(bquote('Number of fixations')) +
  xlab("Time (s)") +  
  #xlab("") +
  scale_x_continuous(limits = c(0,10),breaks = seq(0,10,2)) +
  scale_y_continuous(limits = c(2.0,3.5),breaks = seq(2.0,3.5,0.5)) +
  ggtitle("Fixation number") +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        title=element_text(size=14))

y1 <- ggplot_build(fixnum_plot)$layout$panel_params[[1]]$y$limits[2]
y2 <- y1 - diff(ggplot_build(fixnum_plot)$layout$panel_params[[1]]$y$limits)/25
y3 <- y2 - diff(ggplot_build(fixnum_plot)$layout$panel_params[[1]]$y$limits)/25
y4 <- y3 - diff(ggplot_build(fixnum_plot)$layout$panel_params[[1]]$y$limits)/25

# Add results of post-hoc tests
fixnum_plot <- add_cluster_analysis.actpass(fixnum_plot, fnum_clusters_within, "fnum.high.act_vs_pass", "#FF0000", 1)
fixnum_plot <- add_cluster_analysis.actpass(fixnum_plot, fnum_clusters_within, "fnum.low.act_vs_pass", "#000099", 1)
fixnum_plot <- add_cluster_analysis.highlow(fixnum_plot, fnum_clusters_within, "fnum.passive.high_vs_low", "#b74ac0", 1)
fixnum_plot <- add_cluster_analysis.highlow(fixnum_plot, fnum_clusters_within, "fnum.active.high_vs_low", "#b74ac0", 1)
fixnum_plot <- fixnum_plot +
  annotate("text", x=0, y=y1, label= "HI active vs. HI passive", size = 3, color = "#FF0000", hjust = 0)+
  annotate("text", x=0, y=y2, label= "HI active vs. LI active", size = 3, color = "#b74ac0", hjust = 0)+
  annotate("text", x=0, y=y3, label= "HI passive vs. LI passive", size = 3, color = "#b74ac0", hjust = 0)+
  annotate("text", x=0, y=y4, label= "LI active vs. LI passive", size = 3, color = "#000099", hjust = 0)

fixnum_plot

#ggsave(file = paste(path,"fixnum.png",sep=""), dpi = 600, width = 9, height = 5, units = "in", bg="white")


# Plot Duration of Fixations

fixdur_plot <- fdur %>%
  group_by(bin,cue) %>%
  summarise(mean = mean(value, na.rm= TRUE), se = sd(value,na.rm=TRUE) / sqrt(n())) %>%
  ggplot(aes(x = bin, group = cue, y = mean)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE")+
  geom_line(aes(colour = cue, linetype = cue), size = 0.8) +
  scale_linetype_manual(name="Trial type", values=c( "dotdash","solid",  "dotdash", "solid"), labels = c("HI active","HI passive", "LI active","LI passive"))+
  scale_color_manual(name="Trial type", values = c("#FF0000","#FF0000","#000099", "#000099"),labels = c("HI active","HI passive", "LI active","LI passive")) + 
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = cue), alpha = 0.1) + 
  scale_fill_manual(values = c("#FF0000","#FF0000","#000099", "#000099")) +
  guides(fill=FALSE) +
  ylab(bquote('Fixation duration (ms)')) +
  xlab("Time (s)") +  
  #xlab("") +
  scale_x_continuous(limits = c(0,10),breaks = seq(0,10,2)) +
  scale_y_continuous(limits = c(250,1100),breaks = seq(300,1100,150)) +
  ggtitle("Fixation duration") +
  theme_classic() +
  theme(legend.position = "none")

y1 <- ggplot_build(fixdur_plot)$layout$panel_params[[1]]$y$limits[2]
y2 <- y1 - diff(ggplot_build(fixdur_plot)$layout$panel_params[[1]]$y$limits)/25
y3 <- y2 - diff(ggplot_build(fixdur_plot)$layout$panel_params[[1]]$y$limits)/25
y4 <- y3 - diff(ggplot_build(fixdur_plot)$layout$panel_params[[1]]$y$limits)/25

# Add results of post-hoc tests
fixdur_plot <- add_cluster_analysis.actpass(fixdur_plot, fdur_clusters_within, "fdur.high.act_vs_pass", "#FF0000", 1)
fixdur_plot <- add_cluster_analysis.actpass(fixdur_plot, fdur_clusters_within, "fdur.low.act_vs_pass", "#000099", 1)
fixdur_plot <- add_cluster_analysis.highlow(fixdur_plot, fdur_clusters_within, "fdur.passive.high_vs_low", "#b74ac0", 1)
fixdur_plot <- add_cluster_analysis.highlow(fixdur_plot, fdur_clusters_within, "fdur.active.high_vs_low", "#b74ac0", 1)
# fixdur_plot <- fixdur_plot +
#   annotate("text", x=0, y=y1, label= "HI active vs. HI passive", size = 3, color = "#FF0000", hjust = 0)+
#   annotate("text", x=0, y=y2, label= "HI active vs. LI active", size = 3, color = "#b74ac0", hjust = 0)+
#   annotate("text", x=0, y=y3, label= "HI passive vs. LI passive", size = 3, color = "#b74ac0", hjust = 0)+
#   annotate("text", x=0, y=y4, label= "LI active vs. LI passive", size = 3, color = "#000099", hjust = 0)


fixdur_plot


# Plot Center Bias

cbias_plot <- cb %>%
  group_by(bin,cue) %>%
  summarise(mean = mean(value, na.rm= TRUE), se = sd(value,na.rm=TRUE) / sqrt(n())) %>%
  ggplot(aes(x = bin, group = cue, y = mean)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE")+
  geom_line(aes(colour = cue, linetype = cue), size = 0.8) +
  scale_linetype_manual(name="Trial type", values=c( "dotdash","solid",  "dotdash", "solid"), labels = c("HI active","HI passive", "LI active","LI passive"))+
  scale_color_manual(name="Trial type", values = c("#FF0000","#FF0000","#000099", "#000099"),labels = c("HI active","HI passive", "LI active","LI passive")) + 
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = cue), alpha = 0.1) + 
  scale_fill_manual(values = c("#FF0000","#FF0000","#000099", "#000099")) +
  guides(fill=FALSE) +
  ylab(bquote('Distance from center (px)')) +
  #xlab("") +
  xlab("Time (s)") +
  scale_x_continuous(limits = c(0,10),breaks = seq(0,10,2)) +
  scale_y_continuous(limits = c(75,180),breaks = seq(70,180,20)) +
  ggtitle("Central Bias") +
  theme_classic() +
  theme(legend.position = "none")

y1 <- ggplot_build(cbias_plot)$layout$panel_params[[1]]$y$limits[2]
y2 <- y1 - diff(ggplot_build(cbias_plot)$layout$panel_params[[1]]$y$limits)/25
y3 <- y2 - diff(ggplot_build(cbias_plot)$layout$panel_params[[1]]$y$limits)/25
y4 <- y3 - diff(ggplot_build(cbias_plot)$layout$panel_params[[1]]$y$limits)/25

# Add results of post-hoc tests

cbias_plot <- add_cluster_analysis.actpass(cbias_plot, cb_clusters_within, "cb.high.act_vs_pass", "#FF0000", 1)
cbias_plot <- add_cluster_analysis.actpass(cbias_plot, cb_clusters_within, "cb.low.act_vs_pass", "#000099", 1)
cbias_plot <- add_cluster_analysis.highlow(cbias_plot, cb_clusters_within, "cb.passive.high_vs_low", "#b74ac0", 1)
cbias_plot <- add_cluster_analysis.highlow(cbias_plot, cb_clusters_within, "cb.active.high_vs_low", "#b74ac0", 1)
# cbias_plot <- cbias_plot +
#   annotate("text", x=0, y=y1, label= "HI active vs. HI passive", size = 3, color = "#FF0000", hjust = 0)+
#   annotate("text", x=0, y=y2, label= "HI active vs. LI active", size = 3, color = "#b74ac0", hjust = 0)+
#   annotate("text", x=0, y=y3, label= "HI passive vs. LI passive", size = 3, color = "#b74ac0", hjust = 0)+
#   annotate("text", x=0, y=y4, label= "LI active vs. LI passive", size = 3, color = "#000099", hjust = 0)

cbias_plot



# Plot Alpha

alpha_long <- read.csv2("alpha_long.csv")%>%
  group_by(ID, time, ActPass, Intensity)%>%
  summarize(
    value = mean(value)
  )%>%
  mutate(
    ID = as.factor(ID),
    cue = as.factor(case_when(ActPass == "Passive" & Intensity == "High" ~ 1,
                    ActPass == "Active" & Intensity == "High" ~ 2,
                    ActPass == "Passive" & Intensity == "Low" ~ 3,
                    ActPass == "Active" & Intensity == "Low" ~ 4)),
    time = as.numeric(time)/2-0.25)

alpha_long$cue <- factor(alpha_long$cue, levels = c("1", "2", "3", "4"))

    
alpha_plot <- alpha_long %>%
  group_by(time,cue) %>%
  summarise(mean = mean(value, na.rm= TRUE), se = sd(value,na.rm=TRUE) / sqrt(n())) %>%
  ggplot(aes(x = time, group = cue, y = mean)) +
  geom_rect(ymax=Inf,ymin=-Inf,xmax=10,xmin=2,fill="#EEEEEE")+
  geom_line(aes(colour = cue, linetype = cue), size = 0.8) +
  scale_linetype_manual(name="Trial type", values=c("solid", "dotdash", "solid", "dotdash"), labels = c("HI passive","HI active", "LI passive","LI active"))+
  scale_color_manual(name="Trial type", values = c("#FF0000","#FF0000","#000099", "#000099"),labels = c("HI passive","HI active", "LI passive","LI active")) + 
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se, fill = cue), alpha = 0.1) + 
  scale_fill_manual(values = c("#FF0000","#FF0000","#000099", "#000099")) +
  guides(fill=FALSE) +
  ylab(bquote('Signal Change in %')) +
  #xlab("") +
  xlab("Time (s)") +
  scale_x_continuous(limits = c(-1,10),breaks = seq(0,10,2)) +
  scale_y_continuous(limits = c(-35,10)) +
  ggtitle("Alpha Power Suppression") +
  theme_classic()+
  theme(legend.position = "none")

y1 <- ggplot_build(alpha_plot)$layout$panel_params[[1]]$y$limits[2]
y2 <- y1 - diff(ggplot_build(alpha_plot)$layout$panel_params[[1]]$y$limits)/25
y3 <- y2 - diff(ggplot_build(alpha_plot)$layout$panel_params[[1]]$y$limits)/25
y4 <- y3 - diff(ggplot_build(alpha_plot)$layout$panel_params[[1]]$y$limits)/25

# Add results of post-hoc tests

alpha_plot <- add_cluster_analysis.actpass(alpha_plot, alpha_clusters_within, "alpha.high.act_vs_pass", "#FF0000", 0.5)
alpha_plot <- add_cluster_analysis.actpass(alpha_plot, alpha_clusters_within, "alpha.low.act_vs_pass", "#000099", 0.5)
alpha_plot <- add_cluster_analysis.highlow(alpha_plot, alpha_clusters_within, "alpha.passive.high_vs_low", "#b74ac0", 0.5)
alpha_plot <- add_cluster_analysis.highlow(alpha_plot, alpha_clusters_within, "alpha.active.high_vs_low", "#b74ac0", 0.5)
# alpha_plot <- alpha_plot +
#   annotate("text", x=-1, y=y1, label= "HI active vs. HI passive", size = 3, color = "#FF0000", hjust = 0)+
#   annotate("text", x=-1, y=y2, label= "HI active vs. LI active", size = 3, color = "#b74ac0", hjust = 0)+
#   annotate("text", x=-1, y=y3, label= "HI passive vs. LI passive", size = 3, color = "#b74ac0", hjust = 0)+
#   annotate("text", x=-1, y=y4, label= "LI active vs. LI passive", size = 3, color = "#000099", hjust = 0)
alpha_plot


empty_plot <- ggplot() + theme_void()


all_eyes <- fixdur_plot + cbias_plot/empty_plot 
all_eyes


all_plot <- (hr_plot + sc_plot + pupil_plot) / (cbias_plot + fixdur_plot)
all_plot

all_plot <- cbias_plot + fixdur_plot + 
  alpha_plot + pupil_plot + 
  sc_plot + hr_plot +
  plot_layout(ncol = 2)+
  plot_annotation(tag_levels = 'A')

all_plot
  
#ggsave(file = paste(path,"allplots.png",sep=""), dpi = 600, width = 10, height = 12, units = "in", bg="white")
#ggsave(file = paste(savepath,"allplots.pdf",sep=""), dpi = 600, width = 4, height = 12, units = "in", bg="white")

