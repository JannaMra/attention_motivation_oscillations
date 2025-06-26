############################################################################
# ShockScanning_Flight Oscillations project
# 
# - Score baselines and determine baseline quality
#
library(ggplot2)


path <- paste0(getwd(),"/")
path.prot <- "Data/prot/" %>% paste0(path, .)
savepath <- "Data/Tobii/BL/" %>% paste0(path, .)

exclusions <- c("vp05","vp22", "vp30","vp51") #vp05 & vp30 responded too early; vp51 underwent block 1 twice, vp 22 too many responses in passive trials (empty cells)


# Baseline between -300 and 0 ms relative to stimulus onset
blst <- -300; blen <- 0
# Stimulus duration
dur <- 10000

outlierLimit.eye <- 0.5

# Plot individual baselines to jpg?
plotBL1 <- TRUE
plotBL2 <- TRUE

# Daten laden
msg <- read.table(paste(path,"Data/Tobii/Messages.txt",sep=""))
names(msg) <- c("vp","trial","time")
fixa <- read.table(paste(path,"Data/Tobii/Fixations.txt", sep=""))

# Determine which subjects should be analyzed
vpn = fixa$subject %>% unique() %>% sort() %>% as.character() #all subjects in fixations
vpn.n = length(vpn)
vpn = vpn[vpn %in% exclusions == F] #minus a priori exclusions
vpn.n = length(vpn)

# Exclusions:
vpn <-  vpn[!(vpn %in% c())]


# Iterative outlier removal
outlier_remove <- function(x,sdmult=3) {
  # Remove outlier iteratively
  insample <- rep(1,length(x)); insample[is.na(x)] <- 0
  ok <- FALSE
  while (!ok) {
    xminpos <- (1:length(x))[x==min(x[insample==1])]
    xminpos <- xminpos[xminpos %in% (1:length(insample))[insample==1]][1]
    xmaxpos <- (1:length(x))[x==max(x[insample==1])]
    xmaxpos <- xmaxpos[xmaxpos %in% (1:length(insample))[insample==1]][1]
    tempinsample <- insample; tempinsample[c(xminpos,xmaxpos)] <- 0
    subx <- x[tempinsample==1]
    if (x[xminpos]<(mean(subx)-sdmult*sd(subx))) {
      insample[xminpos] <- 0
      out1 <- TRUE
    } else {
      out1 <- FALSE
    }
    
    if (x[xmaxpos]>(mean(subx)+sdmult*sd(subx))) {
      insample[xmaxpos] <- 0
      out2 <- TRUE
    } else {
      out2 <- FALSE
    }
    
    if (!out1 & !out2) { ok <- TRUE }
  }
  return(insample)
}


bl_quality <- data.frame()
validfix_quality <- data.frame()

# Loop over subjects
for (vp in vpn) {
  #vp= "vp10"
  
  vpcode <- vp #paste0("vp",ifelse(as.numeric(vpn)<10, paste(c("0", vpn), collapse = ""),vpn))
  prot <- read.csv2(paste(path.prot,vpcode,".csv",sep=""),header=TRUE)
  prot$blx  <- NA
  prot$bly  <- NA
  prot$blok1 <- NA # max. 3 SD
  prot$blok2 <- NA # max. spread in px
  prot$blok <- NA  # both criteria combined 
  
  blxok_1 <- rep(0,80)
  blyok_1 <- rep(0,80)
  blxok_2 <- rep(0,80)
  blyok_2 <- rep(0,80)
  
  baseline <- numeric()  # Position (x/y)
  validfix <- numeric()  # Valid fixation time (ms), excluding 1st fixation
  code <- vp
  print(code)
  
  ###########################################
  # 1. Determine baselines
  # Generate empty data field to store data
  
  # Determine trial number
  vpfix <- fixa[tolower(fixa$subject)==code,]
  #ntrial <- nrow(prot)
  
  # Loop over trials to determine trial-by-trial baselines
  for(blockid in 1:2){
    print(blockid)
    baseline <- numeric()
    if(blockid == 1) {
      trialid = 1
      ntrial = 80
    }
    if(blockid == 2) {
      trialid = 81
      ntrial = 160
    }
       
    for(trialid in trialid:ntrial){
    # Select trial data
      fixblock <- fixa%>%
        filter(trial == trialid &
               subject == code &
               block == blockid)
      msgblock <- msg [
        tolower(msg$vp)==code &
          msg$trial==trialid,]
    
    # Check if correct Stimulus marker is present in MSG file
    # if (msgblock$event!=paste(prot$pic[trial],
    #                           ".jpg",sep="")) {
    #   print(paste("Stimulus error: Trial:",trial," Event:",msgblock$event,sep=""))
    # }
    
    # Determine onset (in ms)
    onset <- as.numeric(msgblock$time)
    
    # Subtract onset from timestamps
    fixblock$start  <- as.numeric(fixblock$start)-onset
    fixblock$end <- as.numeric(fixblock$end)-onset
    
    # Cut last fixation to presentation time
    fixblock$end[nrow(fixblock)] <- 
      ifelse(
        fixblock$end[
          nrow(fixblock)]>dur,
        dur,
        fixblock$end[
          nrow(fixblock)])
    
    
    # Calculate valid fixation time
    fixblock$dur <- as.numeric(fixblock$end) - as.numeric(fixblock$start)
    validfix <- c(validfix,sum(fixblock$dur[fixblock$start>0],na.rm=TRUE))
    
    # Calculate baseline as weighted average of fixations
    fixblockbl <- fixblock[fixblock$end>blst & fixblock$start<blen,]
    if (nrow(fixblockbl)>0) {
      # Restrict fixation data to baseline
      fixblockbl$start[1] <- ifelse(head(fixblockbl$start,1)<blst,blst,head(fixblockbl$start,1))
      fixblockbl$end[nrow(fixblockbl)] <- ifelse(tail(fixblockbl$end,1)>blen,blen,tail(fixblockbl$end,1))
      
      
      # Calculate baseline coordinates
      xbl <- sum(fixblockbl$x*fixblockbl$dur)/sum(fixblockbl$dur)
      ybl <- sum(fixblockbl$y*fixblockbl$dur)/sum(fixblockbl$dur)
      
      # Store values
      baseline <- rbind(baseline,c(xbl,ybl))
    } else {
      # When no valid fixations are available store NA as baseline for current trial
      baseline <- rbind(baseline,c(NA,NA))
    }
    }
    if(blockid == 1) {
      baseline_1 <- baseline}
    if(blockid == 2) {
      baseline_2 <- baseline}
  }
  
  # Determine outlier (only perform baseline validation if less than 100 NAs/ at least 10 valid fixations)
  
  if (!sum(is.na(baseline_1))> 100) {
    blxok_1 <- outlier_remove(baseline_1[,1])
    blyok_1 <- outlier_remove(baseline_1[,2])
  } else {print("Too little values for Baseline Validation in Block 1")}
  
  if (!sum(is.na(baseline_2))> 100) {
    blxok_2 <- outlier_remove(baseline_2[,1])
    blyok_2 <- outlier_remove(baseline_2[,2])
  } else {print("Too little values for Baseline Validation in Block 2")}
  
  # Baseline is valid when x and y coordinates are ok (i.e. no outlier)
  blok_1 <- as.numeric((blxok_1==1) & (blyok_1==1))
  blok_2 <- as.numeric((blxok_2==1) & (blyok_2==1))
  
  ##NEW AM Error: in $ data frame ... replacement has 183 rows, data has 60
  #prot$blx <- cut(prot$blx, baseline[,1])
  #prot$bly <- baseline[,2]
  
  prot$blx <- c(baseline_1[,1], baseline_2[,1])
  prot$bly <- c(baseline_1[,2], baseline_2[,2])
  prot$blok1 <- c(blok_1, blok_2) # outlier because more than 3 SD from mean 
  prot$validfix <- validfix
  # hier weiter
  meanx1 <- rep(mean(baseline_1[,1],na.rm = TRUE), 80)
  meanx2 <- rep(mean(baseline_2[,1],na.rm = TRUE), 80)
  meanx <- c(meanx1, meanx2)
  prot$blmeanx <- meanx
  meany1 <- rep(mean(baseline_1[,2],na.rm = TRUE), 80)
  meany2 <- rep(mean(baseline_2[,2],na.rm = TRUE), 80)
  meany <- c(meany1, meany2)
  prot$blmeany <- meany

  prot <- prot %>%
    mutate(blok2 = case_when(abs(blmeanx - blx) > 75 | abs(blmeany - bly) > 75 ~ 0, #outlier because more than 150 pixels from mean
                             abs(blmeanx - blx) < 75 & abs(blmeany - bly) < 75 ~ 1
    )
    )
  prot <- prot %>%
    mutate(blok = case_when(blok1 == 0 | blok2 == 0 ~ 0, # outlier because either more than 3 SD or more than 150 pixels away
                            blok1 == 1 & blok2 == 1 ~ 1))
  blok <- prot$blok
  
  prot1 <- head(prot,80)
  prot2 <- tail(prot,80)
  
  # Testplot
  
  if (plotBL1) {
    png(paste(savepath,vpcode,"_1",".jpg",sep=""),width=500,height=500,pointsize=18)
    plot(prot1$blx,prot1$bly,pch=16,col="black",xlab="x (px)",ylab="y (px)",xlim=c(0,1280),ylim=c(0,1024))
    points(prot1%>% filter(blok==0) %>% select(blx,bly),pch=16,col="red")
    title(paste0(code,"_1"))
    dev.off()
  }
  if (plotBL2) {
    png(paste(savepath,vpcode,"_2",".jpg",sep=""),width=500,height=500,pointsize=18)
    plot(prot2$blx,prot2$bly,pch=16,col="black",xlab="x (px)",ylab="y (px)",xlim=c(0,1280),ylim=c(0,1024))
    points(prot2%>% filter(blok==0) %>% select(blx,bly),pch=16,col="red")
    title(paste0(code,"_2"))
    dev.off()
  }
  
  # Store number of valid baselines per subject
  bl_quality_vp <- data.frame(nrow(prot),sum(blok))
  
  # Store number of percent valid fixations during trial per subject 
  validfix_vp <- data.frame(nrow(baseline),mean(validfix)/10000)
  
  # Store spread of baselines in x and y direction
  # xrng <- max(baseline[blok==1,1])-min(baseline[blok==1,1])
  # yrng <- max(baseline[blok==1,2])-min(baseline[blok==1,2])
  # bl_quality_vp <- c(bl_quality_vp,xrng,yrng)
  
  write.csv2(prot,paste(path.prot,vpcode,"_BL.csv",sep=""),row.names=FALSE,quote=FALSE)
  
  bl_quality <- rbind(bl_quality,bl_quality_vp)
  validfix_quality <- rbind(validfix_quality, validfix_vp)
  
}

erg <- data.frame(vpn,bl_quality,validfix_quality[2],row.names=NULL)
names(erg) <- c("code","alltrials","blok","validfix")

# Mark problematic cases
erg$pblok <- erg[,3]/erg[,2]
erg$prob <- ifelse(erg$pblok<= 0.5,1,0)
sum(erg$prob)
#erg$problem  <- as.numeric((erg$xrng>150) | (erg$yrng>150))
write.csv2(erg,paste(path.eye,"Results_BaselineCheck_blockwise.csv",sep=""),row.names=FALSE,quote=FALSE)

# Define cases with too many outliers (>50%) to exclude from further analyses
eye.invalid.bl <- erg %>% filter(prob == 1) %>% pull(code)
#eye.invalid.bl <- c("vp04", "vp08", "vp09", "vp20", "vp31", "vp46")
  
  
#erg %>% group_by(code) %>% with(hist(pblok, breaks=20, main="Valid baselines gesamt")); abline(v = outlierLimit.eye, col="red", lwd=2, lty=2)
plot <- erg %>%
  ggplot(aes(x=pblok)) +
  geom_histogram(binwidth=0.1, color = "#e9ecef", alpha=0.4, position = 'identity') +
  geom_vline(xintercept= 0.5, color = "red") +
  xlab("Percent valid Baselines")

print(plot)


