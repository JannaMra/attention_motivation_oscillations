##################################################################
# Shockscanning Flight Choice
# Project
#
# Average and plot fixation density maps for each condition
#
# Cue-Kodierung: 0 = Flight, 1 = Shock, 2 = Safety

#rm(list=ls())

require(viridis)  # For colormaps
require(png)

#rm(list=ls())

computeAV <- TRUE
plotPNG <- TRUE

zmin <- -2
zmax <- 0

logscale <- TRUE

plot.heatmap <- function(heatm,...) {
  # transpose and flip matrix for printing
  plotmap <- t(heatm)
  plotmap <- plotmap[,ncol(plotmap):1]
  
  image(1:nrow(plotmap),1:ncol(plotmap),plotmap,
        xlim=c(0,nrow(plotmap)),ylim=c(0,ncol(plotmap)),asp=1,axes=FALSE,...)
  
  #axis(1,seq(0,nrow(plotmap),200),seq(0,nrow(plotmap),200)-nrow(plotmap)/2,pos=0)
  #axis(2,seq(0,ncol(plotmap),200),(-1)*(seq(0,ncol(plotmap),200)-ncol(plotmap)/2),las=2,pos=0)
}

#jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
jet.colors <- colorRampPalette(c("#00007F", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

#path <- "C:/Users/jat41gk/Documents/Projekte/Threat Oscillations/Threat Oscillations/"

protpath <- "Data/prot/" %>% paste0(path, .)
etpath <- "Data/Plots/Data/" %>% paste0(path, .) 
savepath <- "Data/Plots/" %>% paste0(path, .)

# Exclusions:
#print(eye.invalid.bl)

# Determine which subjects should be analyzed
vpn <- list.files(etpath) 
#vpn <- vpn[grepl("^tfo", vpn)]

vps=vpn

if (computeAV) {
  for (vpn in vps) {
    #vpn <- "1"
    print(vpn)
    #vpn <- as.numeric(vpn)
    #vpcode <- paste0("tfo",ifelse(vpn<10, paste(c("0", vpn), collapse = ""),vpn))    #ifelse(vp<10,0,""),vp)
    #print(vpcode)
    
    load(paste(etpath,vpn,sep=""))
    
    if (vpn==vps[1]) {
      fixmaps <- fixdens
    } else {
      for (i in 1:length(fixdens)) {
        fixmaps[[i]] <- fixmaps[[i]]+fixdens[[i]]
      }
    }
  }
  
  for (i in 1:length(fixmaps)) {
    fixmaps[[i]] <- fixmaps[[i]]/length(vps)
  }
  
  save(fixmaps, file=paste(etpath,"All_subjects.RData",sep=""))
} else {
  load(file=paste(etpath,"All_subjects.RData",sep=""))
}

# Plot density maps
if (plotPNG) {
  png(filename=paste(savepath,"FixationDensity_Cues.png",sep=""), width=4000, height=1500, res=300, pointsize=14)  
} else {
  x11()  
}

layout(matrix(c(1,2,3,4,5),nrow=1,byrow=TRUE),widths=c(5,5,5,5,2))

# Plot heatmaps and determine value range
par(mar=c(1,1,1,1))
zrange <- numeric()
for (i in 1:4) {
  plotmap <- fixmaps[[i]]
  zrange <- rbind(zrange, range(range(fixmaps[[i]])))
  
  if (logscale) {
    # logarithmic
    colorBreaks = c(0,10^(seq(zmin, zmax, (zmax-zmin)/100))) # set color breaks however you want
  } else {  
    # linear
    colorBreaks = seq(0, 10^zmax, (10^zmax)/100) # set color breaks however you want
  }
  
  plot.heatmap(plotmap,breaks=colorBreaks,col=plasma(length(colorBreaks)-1),xlab="",ylab="",zlim=c(0,1))
  #plot.heatmap(plotmap,col=plasma(255),xlab="",ylab="",zlim=c(0,1))
  #abline(v=ncol(plotmap)/2,h=nrow(plotmap)/2,col="white",lty="dashed")
}

# Colorbar
zlimval <- c(zmin,zmax)

par(mar=c(0,0,0,0))
plot(1,1,type="n",xlim=c(-2,3),ylim=c(-1,2),xlab="",ylab="",axes=FALSE)
lut <- plasma(255)
scale = (length(lut)-1)
for (j in 1:(length(lut)-1)) {
  y = (j-1)/scale
  rect(0,y,0.7,y+1/scale, col=lut[j], border=NA)
}

# Axis
axis(2,c(0,1),as.character(c(0,10^zmax)),pos=0,lwd=0,las=2)
rect(0,0,0.7,1,col="NA",border="black",lwd=1)

# Ticks
if (logscale) {
  # logarithmic
  colorBreaks = c(0,10^(seq(zmin, zmax, (zmax-zmin)/3)))
} else {
  # linear
  colorBreaks = seq(0, 10^zmax, (10^zmax)/5)
}

for (j in 1:length(colorBreaks)) {
  segments(-0.3,colorBreaks[j],0,colorBreaks[j],col="black")
}

if (plotPNG) {
  dev.off()
}

