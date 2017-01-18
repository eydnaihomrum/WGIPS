PlotNascMapBubbles <- function(StoXbaseline,spp,SAres = 5) {

### function to create NASC bubble plot from StoX baseline report
###
### input: -StoXbaseline: directory where baseline report is stored
###        -spp: 3 capital letter species code with apostrophes [e.g. "HER"]
###        -SAres: resolution of plot in nmi

### load required packages
    library(mapdata)

### get interval data summed over depth channels
    setwd(StoXbaseline)
    spp_data_int <- read.table("4_NASC_NASC.txt",header=TRUE)
    spp_data_int <- spp_data_int[c("SampleUnit","NASC")]
    log <- unlist(strsplit(as.character(spp_data_int$EDSU),"/"))
    logpos <- seq(2,length(log),4)
    logs <- as.numeric(log[logpos])
    nasc <- spp_data_int$NASC
    nascinfo <- as.data.frame(cbind(logs,nasc))
    names(nascinfo) <- c("LOG","SA")
    loginfo <- read.table("3_FilterAcoustic_AcousticData_DistanceFrequency.txt",sep="\t",header=TRUE)
    loginfo <- loginfo[c("log_start","lon_start","lat_start")]
    names(loginfo) <- c("LOG","ACLON","ACLAT")
    loginfo <- loginfo[!duplicated(loginfo),]
    spp_data_int <- merge(nascinfo,loginfo,all=TRUE)
    spp_data_int$SA[is.na(spp_data_int$SA)] <- 0
    spp_data_int <- spp_data_int[order(spp_data_int$LOG),]

### sub-sample by SAres interval
    int <- SAres
    int_steps <- seq(ceiling(int/2),length(spp_data_int$LOG),int)
    int_data <- matrix(NA,length(int_steps),3)
    for(step in 1:length(int_steps)){
      int_data[step,1] <- spp_data_int$ACLON[int_steps[step]]
      int_data[step,2] <- spp_data_int$ACLAT[int_steps[step]]
      int_data[step,3] <- mean(spp_data_int$SA[(step*int-(int-1)):(step*int)])
    }
    
### bubble plot
    minlon <- min(int_data[,1],na.rm=TRUE)
    maxlon <- max(int_data[,1],na.rm=TRUE)
    minlat <- min(int_data[,2],na.rm=TRUE)
    maxlat <- max(int_data[,2],na.rm=TRUE)
    extfact <- 0.3
    x1 <- round((minlon - extfact*(maxlon-minlon))/1)*1
    x2 <- round((maxlon + extfact*(maxlon-minlon))/1)*1
    y1 <- round((minlat - extfact*(maxlat-minlat))/0.5)*0.5
    y2 <- round((maxlat + extfact*(maxlat-minlat))/0.5)*0.5
    
    tiff(file=paste(spp,"_SAplot.tiff",sep=""),width=10,height=(y2-y1)*2/(x2-x1)*10,units="cm",res=300)
    par(mfrow=c(1,1), mar=c(4,4,2,0.5), oma=c(1,1,1,1),family="sans", font=2, xaxs="i",yaxs="i",bty="n")
    
    maxSA <- max(int_data[,3],na.rm=TRUE)
    plot(int_data[,1],
         int_data[,2],
         col="white",
         xlim=c(x1,x2),
         ylim=c(y1,y2),
         yaxs="i",
         xaxs="i",
         xlab="",
         ylab="",
         axes=FALSE)
    grid(nx=(x2-x1)/1,ny=(y2-y1)/0.5,lty=1)
    symbols(int_data[,1],
            int_data[,2],
            circles=sqrt(int_data[,3]/maxSA),
            inches=0.1,
            lwd=1,
            yaxs="i",
            xaxs="i",
            xlim=c(x1,x2),
            ylim=c(y1,y2),
            xlab="",
            ylab="",add=TRUE)

    map('worldHires',add=TRUE,col=colors()[226],fill=T)
    
    #title(paste(unique(acoustic_data$CRUISE),", NASC for species:  ",spp,sep=""),cex.main=0.75)
    box(bty="O",lwd=1) # bty=box type to surround plot (e.g. "O", "L", "7", "C", "U" to set box shape or "n" for none)
    axis(side = 1, at = c(seq(x1,x2,1)), labels = as.character(format(seq(x1,x2,1), nsmall = 0)),las=1, cex.axis = 0.5,lwd=1,font=1,mgp=c(3,0.5,0))
    mtext("Longitude", side=1,line=2, cex=1,las=1)
    axis(side = 2, at = seq(y1,y2,0.5), labels = as.character(format(seq(y1,y2,0.5), nsmall = 1)),las=1, cex.axis = 0.5,lwd=1,font=1,mgp=c(3,1,0))
    mtext("Latitude", side=2,line=3, cex=1,las=3)
    
    dev.off()     

}

