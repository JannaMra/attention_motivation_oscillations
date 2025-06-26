##################################################################
# Shockscanning Flight Choice
# Projekt 
#
# Calculate fixation density maps for each condition and observer
#
# Design:
# cue 1: high intensity shock (shock)
# cue 2: high intensity flight (flight)
# cue 3: low intensity shock (safety)
# cue 4: low intensity flight (active control)

#rm(list=ls())

require("spatstat")
#require("readbitmap")
#require("jpeg")
#require("png")

path <- paste0(getwd(),"/")
protpath <- "Data/prot/" %>% paste0(path, .)
savepath <- "Data/Plots/Data/" %>% paste0(path, .)


st <- 0; en <- 8000  # Scoring range

# Baseline between -300 and 0 ms relative to stimulus onset
blst <- -300; blen <- 0

# Image dimensions:
# screen resolution: 1920 x 1200 pixels
screen.height = 1024 #height screen in pix
screen.width  = 1280 # width of screen in pix

wsx <- 1280 ; wsy <- 1024 
# Image size: 768 x 576 pixels (visual angle of 24.07 x 16.3)


# Wie viele Pixel result in a visual angle of 1?
# Distance from monitor: 500 mm
# Monitor: BOLDscreen 32 LCD for fMRI
# 1920 x 1080, RGB colour with fixed 120Hz frame rate
# Display size 698.4mm x 392.9mm 
# Resolution: 1920 x 1080 Pixel
# pic size: 768 × 576 Pixel
visdeg <- 54  

# Eyetrackingdaten laden 
msg <- read.table(paste(path,"Data/Tobii/Messages.txt",sep=""))
names(msg) <- c("vp","trial","time")
fixa <- read.table(paste(path,"Data/Tobii/Fixations.txt", sep=""))

# Exclusions:
print(eye.invalid.bl)

# Determine which subjects should be analyzed
#vpn <- c(paste("tfo",ifelse(1:16<10,"0",""),1:16,sep=""))
vpn = fixa$subject %>% unique() %>% sort() %>% as.character() #all subjects in fixations
vpn.n = length(vpn)
vpn = vpn[vpn %in% exclusions == F] #minus a priori exclusions
vpn <-  vpn[!(vpn %in% eye.invalid.bl)]
vpn.n = length(vpn)


# Loop across participants and save individual fixation density maps for each subject
fixinval <- numeric()
#for (vp in vpn) {
for (vp in vpn) {
  #vp <- "vp01"
  print(vp)
  # read in data
  vpcode <- vp #paste0(ifelse(vpn<10, paste(c("0", vpn), collapse = ""),vpn))    #ifelse(vp<10,0,""),vp)
  print(vpcode)
  
  prot <- read.csv2(paste(protpath,vpcode,"_BL.csv",sep=""))
  # prot <- prot %>% mutate(cue = case_when(
  #   cue == "1" ~ "0",
  #   cue == "0" ~ "1",
  #   cue == "2" ~ "2"
  # ))
  
  fixdens <- list()
  fixinvalvp <- numeric()
  
  # loop through different cue types
  for (cue in c(1,2,3,4)) {
    if (exists("fixall")) { rm("fixall") }
    for (trialnr in 1:nrow(prot)) {  
      #trialnr <- 1
      # only analyse when Baseline OK and cue available
      if ((prot$blok[trialnr]== 1) & prot$cue[trialnr]==cue) {
        fixblock <- fixa[fixa$subject==vp & fixa$trial==trialnr,]
        
        # define onset
        onset  <- msg$time[msg$vp==vp & msg$trial==trialnr]

        if (nrow(fixblock)>0) {
          # shift to onset
          fixblock$start  <- fixblock$start-onset - 2000 # time of fixcross and first half of trial
          fixblock$end <- fixblock$end-onset - 2000 # time of fixcross
          
          # consider relevant fixations only
          fixblock <- fixblock[fixblock$start>st & fixblock$end<=en,] # switched start and end here ?

          if (nrow(fixblock)>0) { 
            # correct fixations for baseline
            fixblock$x <- fixblock$x-prot$blx[trialnr]
            fixblock$y <- fixblock$y-prot$bly[trialnr]
            
            if (!exists("fixall")) {
              fixall <- fixblock
            } else {
              fixall <- rbind(fixall, fixblock)
            }
          }
        } 
      }
    }
    
    fixall$dur <- fixall$end - fixall$start
    
    # Distance from center of display
    fixall$dist <- sqrt(fixall$x^2+fixall$y^2)
    
    # Valid sample (i.e. on screen?)
    fixall$valid <- 0
    fixall$valid[(fixall$x>(-wsx/2)) & (fixall$x<(wsx/2)) &
                 (fixall$y>(-wsy/2)) & (fixall$y<(wsy/2))] <- 1
    
    fixinvalvp <- c(fixinvalvp, sum(fixall$valid==0))
    
    # Build data vectors with 1ms sample rows
    xsmp <- rep(fixall$x,fixall$dur)
    ysmp <- rep(fixall$y,fixall$dur)
    
    # Calculate fixation density map
    X <- suppressWarnings(ppp(xsmp+768/2, (-1)*(ysmp-576/2), c(0,768), c(0,576)))
    #suppressWarnings(plot(X))
    Y <- density(X, sigma=visdeg, dimyx=c(576,768)) # from spatsat: sigma(1SD)=visdeg -> M-1SD & M+1SD = 2*visdeg
    # Normalize image
    Yden <- (Y$v-min(Y$v))/(max(Y$v)-min(Y$v))
    # Flip rows
    Yden <- Yden[nrow(Yden):1,]
    
    fixdens[[cue]] <- Yden
  } # end cue loop
  
  fixinval <- rbind(fixinval,fixinvalvp)
  save(fixdens, file=paste(savepath, vpcode, ".RData", sep=""))
} # end of trial loop
  

