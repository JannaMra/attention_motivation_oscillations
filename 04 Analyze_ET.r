##################################################################
# Shockscanning Threat Flight Oscillations
# Design:
# cue 1: high intensity shock (shock)
# cue 2: high intensity flight (flight)
# cue 3: low intensity shock (safety)
# cue 4: low intensity flight (active control)
#
# Score Eye-Tracking data binwise as weighted averages:
# - Distance of fixations to screen center
# - Fixation durations
# - Fixation numbers
# - Blinks

# 
############################################################################


#rm(list=ls())

options(warn=2)  # Turn warnings into errors / 0=default

st <- 0; en <- 8000  # Scoring range
binsize <- 1000
# Bin timepoints
breaks <- seq(st,en,binsize)


path <- paste0(getwd(),"/")
savepath <- "Data/prot/" %>% paste0(path, .)  

# Eyetrackingdaten laden
msg <- read.table(paste(path,"Data/Tobii/Messages.txt",sep=""))
names(msg) <- c("vp","trial","time")
fixa <- read.table(paste(path,"Data/Tobii/Fixations.txt", sep=""))

# Exclusions:
# Exclusions:
exclusions <- c("vp05","vp22","vp30","vp51")
eye.invalid.bl <- c("vp04", "vp08", "vp09", "vp20", "vp31", "vp46")
print(eye.invalid.bl)

# Determine which subjects should be analyzed
vpn = fixa$subject %>% unique() %>% sort() %>% as.character() #all subjects in fixations
vpn.n = length(vpn)
vpn = vpn[vpn %in% exclusions == F] #minus a priori exclusions
vpn <-  vpn[!(vpn %in% eye.invalid.bl)]
vpn.n = length(vpn)


# Loop over subjects
for (vp in vpn) {
  
  #vp <- "vp01"
  #vpn <- as.numeric(vpn)
  vpcode <- vp
  print(vpcode)
  
  prot <- read.csv2(paste(path.prot,vpcode,"_BL.csv",sep=""))
  
  # Correct baselines of invalid trials to means of valid trials
  prot$blx[prot$blok==0] <- prot$blmeanx[prot$blok==0] #mean(prot$blx[prot$blok==1])
  prot$bly[prot$blok==0] <- prot$blmeany[prot$blok==0] #mean(prot$bly[prot$blok==1])
  
  # create progress bar
  pb <- txtProgressBar(min=0, max=nrow(prot), style=3)

  # Determine trial number
  ntrial <- nrow(prot)
  
  # Loop over trials to determine values trial-by-trial
  all_cb <- numeric(); all_fixdur <- numeric(); all_fixn <- numeric(); #all_blinks <- numeric()
  for (trial in 1:ntrial) {
    #trial <- 1
    # Set progress bar
    setTxtProgressBar(pb, trial)
    
    # Select trial data
    fixblock <- fixa[tolower(fixa$subject)==vp & fixa$trial==trial,]
    #sacblock <- sac[tolower(sac$vp)==vp & sac$trial==trial,]
    msgblock <- msg[tolower(msg$vp)==vp & msg$trial==trial,]
    
    # Check if correct Stimulus marker is present in MSG file
    # if (msgblock$event!=paste("Stimulus ",prot$pic[trial],".png",sep="")) {
    #   print(paste("Stimulus error: Trial:",trial," Event:",msgblock$event,sep=""))
    # }
    
    # Determine onset (in ms)
    onset <- msgblock$time
    
    # Subtract onset from timestamps and remove cue period
    fixblock$start  <- fixblock$start-onset - 2000
    fixblock$end <- fixblock$end-onset - 2000
    #sacblock$start  <- sacblock$start-onset - 2000
    #sacblock$end <- sacblock$end-onset - 2000
    
    # Extract relevant stimulation period (exclude first fixation if it extends beyond the scoring window)
    fixblock <- fixblock[fixblock$start>st & fixblock$start<=en,]
    #sacblock <- sacblock[sacblock$start>st & sacblock$start<=en,]
    
    # Center bias, fixation duration, fixation number
    if (nrow(fixblock)>0) {
      # Cut last fixation to scoring window
      fixblock$end[nrow(fixblock)] <- ifelse(fixblock$end[nrow(fixblock)]>en, en, fixblock$end[nrow(fixblock)])
      
      # Fixation baseline correction
      fixblock$x <- fixblock$x-prot$blx[trial]
      fixblock$y <- fixblock$y-prot$bly[trial]
      
      # Distance to screen center
      fixblock$dist <- sqrt(fixblock$x^2+fixblock$y^2)
      
      # Fixation durations
      fixblock$duration <- fixblock$end - fixblock$start
      
      # Fixation number
      fixblock$fixn <- 1:nrow(fixblock)
      
      # Valid sample (i.e. on screen?)
      fixblock$valid <- 0
      fixblock$valid[(fixblock$x>=(-screen.width/2)) & (fixblock$x<=(screen.width/2)) & 
                     (fixblock$y>=(-screen.height/2)) & (fixblock$y<=(screen.height/2))] <- 1

      ii <- 0 
      while (ii <= nrow(fixblock)-1){ 
        ii <- ii + 1 
        
        for (bb in 1:length(breaks)){
          
          if((fixblock$start[ii] <= breaks[bb]) && (fixblock$end[ii] > breaks[bb])) {
            fixblock <- rbind(fixblock[1:ii,],fixblock[ii,],fixblock[-(1:ii),]) # replicate row right below previous row
            fixblock$end[ii+1] <- fixblock$end[ii]
            fixblock$end[ii] <- breaks[bb]
            fixblock$start[ii+1] <- breaks[bb]+1
            ii+1
          }
        }
      }
      
      # Calculate duration of new fixation list 
      fixblock$duration.f <- fixblock$end-fixblock$start
      
      # Calculate fractions within each bin
      fixblock$f <- fixblock$duration.f/fixblock$duration
      
      # Add bin numbers
      fixblock$bin <- 1
      for (i in 1:(length(breaks)-1)) {
        fixblock$bin[(fixblock$start>breaks[i]) & (fixblock$end<=breaks[i+1])] <- i
      }
      fixblock$bin <- factor(fixblock$bin, levels=1:(length(breaks)-1))
      
      # Exclude invalid fixations and remove fixations with a duration of 0
      fixblock <- fixblock[(fixblock$valid==1) & (fixblock$duration.f>0),]
      
      if (nrow(fixblock)>0) {
        notrialdata <- FALSE
        
        #######################################
        # Fixation number
        trial_fixn <- tapply(fixblock$f, fixblock$bin, sum)
        trial_fixn[is.na(trial_fixn)] <- 0   # Change NA to 0 fixations available within the bin
        all_fixn <- rbind(all_fixn, trial_fixn)
        
        #######################################
        # Fixation duration
        trial_fixdur <- tapply(fixblock$duration.f, fixblock$bin, sum)/tapply(fixblock$f, fixblock$bin, sum)
        all_fixdur <- rbind(all_fixdur, trial_fixdur)
        
        #######################################
        # Center bias
        trial_cb <- tapply(fixblock$dist*fixblock$duration.f, fixblock$bin, sum)/tapply(fixblock$duration.f, fixblock$bin, sum)
        all_cb <- rbind(all_cb, trial_cb)
      } else {
        notrialdata <- TRUE
      }
    } else {
      notrialdata <- TRUE
    }
    
    if (notrialdata) {
      # When no data are available
      all_fixn   <- rbind(all_fixn, rep(NA,length(breaks)-1))
      all_fixdur <- rbind(all_fixdur, rep(NA,length(breaks)-1))
      all_cb     <- rbind(all_cb, rep(NA,length(breaks)-1))
    }
    
    #######################################
    # Blinks
    # if (nrow(sacblock)>0) {
    #   blinktl <- rep(0,en)
    #   
    #   for (i in 1:nrow(sacblock)) {
    #     if (sacblock$blink[i]=="true") { 
    #       blinktl[sacblock$start[i]:sacblock$end[i]] <- 1
    #     }
    #   }
    #   
    #   # Real time scaling
    #   blinktrial <- numeric()
    #   for (sti in seq(st,en-1,binsize)) {
    #     blinktrial <- c(blinktrial,sum(blinktl[sti:(sti+binsize-1)])/binsize)
    #   }
    #   
    #   all_blinks <- rbind(all_blinks,blinktrial)
    # } else {
    #   # When no data are available
    #   all_blinks <- rbind(all_blinks, rep(NA,length(breaks)-1))
    # }
  }
  
  close(pb)
  
  # Prepare data.frame for saving data
  all_fixn   <- data.frame(all_fixn)
  names(all_fixn) <- paste0(1:(length(breaks)-1),"_FixN")
  all_fixdur <- data.frame(all_fixdur)
  names(all_fixdur) <- paste0(1:(length(breaks)-1),"_FixDur")
  all_cb     <- data.frame(all_cb)
  names(all_cb) <- paste0(1:(length(breaks)-1),"_CB")
  #all_blinks <- data.frame(all_blinks)
  #names(all_blinks) <- paste0(1:(length(breaks)-1),"_Blinks")
  
  out.df <- cbind(prot, all_cb, all_fixn, all_fixdur) #all_blinks)
  
  # Save results (for later aggregation)
  write.csv2(out.df,paste(path.prot,vpcode,"_et.csv",sep=""),row.names=FALSE,quote=FALSE)
}








