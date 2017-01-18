#-------------------------------------------------------------------------------
# Code to estimate uncertainty of transect design outputs from the
# survey design software 'Distance 6.2' (http://distancesampling.org/)
#
# download:
# http://distancesampling.org/Distance/old-versions/downloads/d62setup.exe
#
# By: Sascha Fässler
# IMARES
# 20 April 2016
#-------------------------------------------------------------------------------

rm(list = ls()); graphics.off(); start.time <- Sys.time()
library(RGeostats)
library(PBSmapping)
library(fields)
library(Hmisc)

path <-
  "~/Dropbox/survey demo/survey design"

try(setwd(path),silent = TRUE)

# Load specific geostatistical functions
source("ana.x2y.r")
source("ana.model.r")
source("hermite.r")
source("fn.of.Yp.r")
surveydata.dir <- file.path(".","survey data")
transects.dir <- file.path(".","transects")
polygons.dir <- file.path(".","polygons")
uncertainty.dir <- file.path(".","uncertainty")
europa <- read.table("europa.txt",header = T)

for(str in 1:length(c(91))) {
  
  stratum    <- c(91)[str]

# define variable and fixed input parameters for analysis
EDSU.length   <-    1        #length of the intervals the planned survey [nmi]
transpacing   <-    5        #fixed transect spacing of waypoint input data file in nmi (fixed for now!)
# stratum       <-    91       #stratum numbers
yrs           <-    c(2005,2008:2015)     #survey data year
nbtransimu    <-    250      #define number of transect design simulations
transpacmulti <-    c(1,2,3,4,6) #transect spacing (multitude) scenarioes to be analysed (1 = 5nmi, 2 = 10nmi,...)

yrRSE <-
  array(dim = c(length(1:nbtransimu),3,length(transpacmulti),length(yrs)))
dimnames(yrRSE)[[1]] <- 1:nbtransimu
dimnames(yrRSE)[[2]] <-
  c("simu.mean.sa","sample.mean.sa","sample.RSE.sa")
dimnames(yrRSE)[[3]] <- paste(as.character(transpacing * transpacmulti),"nmi",sep="")
dimnames(yrRSE)[[4]] <- yrs

for (y in 1:length(yrs)) {
  year <- yrs[y]
  # get 5nmi spacing transect design waypoint data (derived from 'Distance software') for the selected stratum
  txectdata <-
    read.table(
      file = file.path(
        transects.dir,paste("stratum ",stratum," 5nm_waypoints.txt",sep = "")
      ),sep = "\t",col.names = c(
        "lon_m","lat_m","lon_dec","lat_dec","transect_ID","stratum_ID"
      )
    )
  strata.select <- as.character(stratum) #select stratum
  polyfile.latlon <-
    list.files(polygons.dir,pattern = paste(paste("_",strata.select,".LatLon",sep =
                                                    ""), collapse = "|"))
  polyfile <-
    read.table(
      file = file.path(polygons.dir,polyfile.latlon),sep = "\t",head = F,col.names = c("lon","lat")
    )
  
  # load historical HERAS data and extract samples inside the polygon
  load(file.path(surveydata.dir,"HERAS.RData"))
  dat.r3 <- data.frame(lon = polyfile$lon,lat = polyfile$lat)
  round.3  <-
    data.frame(
      PID = rep(1,(nrow(dat.r3) + 1)), POS = 1:(nrow(dat.r3) + 1),X = c(dat.r3$lon,dat.r3$lon[1]),Y =
        c(dat.r3$lat,dat.r3$lat[1])
    )
  round.3  <- as.PolySet(round.3,projection = 'UTM')
  events <-
    data.frame(
      EID = 1:nrow(HERAS.data),X = HERAS.data$ACLON,Y = HERAS.data$ACLAT
    ) # Must first creat an events data frame
  events <- as.EventData(events,projection = "UTM")
  r3 <- findPolys(events,round.3,maxRows = 1e+06)
  stratsurv.data <-
    data.frame(
      year = HERAS.data$YEAR[r3$EID],long = HERAS.data$ACLON[r3$EID],lat = HERAS.data$ACLAT[r3$EID],sa =
        HERAS.data$HER[r3$EID]
    )
  
  #conditional geostatistical simulation from the survey year data
  data <- stratsurv.data[which(stratsurv.data$year == year),]
  data <- data[,2:4]
  db.data <- db.create(data,flag.grid = FALSE,ndim = 2,autoname = F)
  pol <-
    read.table(
      file = file.path(
        polygons.dir,paste("stratum_",stratum,".LatLon.txt",sep = "")
      ),sep = "\t",col.names = c("x","y")
    )
  db.poly <- polygon.create(x = pol$x, y = pol$y)
  plot(
    db.data,inches = 5,asp = 1 / cos(mean(data$lat) * pi / 180),title = paste("Acoustic backscatter sampled by the survey",year)
  )
  plot(db.poly,col = 4,add = T)
  polygon(europa,col = 8);box()
  europa.p <- projec.operate(x = europa$x,y = europa$y)
  intperlong <-
    round(30 / EDSU.length) #assuming 30nmi per degree longitude
  res <-
    c(round((max(pol$x) - min(pol$x)) * intperlong),round((max(pol$y) - min(pol$y)) *
                                                            (2 * intperlong))) #grid resolution based on samples of 5nm interval resolution
  db.data <-
    infl(
      db.data,polygon = db.poly,nodes = res,origin = c(min(pol$x),min(pol$y)),extend =
        c(max(pol$x) - min(pol$x),max(pol$y) - min(pol$y)),dmax = 12.5,plot = T,asp =
        1
    ); polygon(europa,col = 8);box()
  # TRANSFORM THE DATA INTO GAUSSIAN: Empirical anamorphosis
  db.data <- ana.x2y(db.data,flag.plot = T)
  db.data <- db.rename(db.data,"gauss","Yp")
  # DEFINING THE MODEL:
  # 1. Modeling the histogram (gaussian anamorphosis i.e. transformation)
  f <- ana.model(db.data,name.Z = "sa",name.Yp = "Yp")
  ycut <- min(db.data[,"Yp"],na.rm = T)
  f.phi <- function(y,ycut = ycut,fct = f) {
    z <- y
    z[y <= ycut] <- 0
    z[y > ycut & !is.na(y)] <- fct(y[y > ycut & !is.na(y)])
    z
  }
  # 2. Modeling the variogram of the gaussian variable
  n.H <- 50
  fn.Yp <-
    fn.of.Yp(ycut,n.H)    ### coefficient of the Hermits polynomial of the function Y --> Yp
  var.Yp.mod <- sum(fn.Yp[-1] ^ 2) ### variance of Yp in the model
  ro <- seq(-0.5,1.22,by = 0.01)
  cov.Yp <- NA
  for (i in 1:length(ro)) {
    cov.Yp[i] <- sum(fn.Yp[-1] ^ 2 * ro[i] ^ seq(1,(n.H - 1)))
  }
  f.Y.vs.Yp <- approxfun(x = cov.Yp,y = ro,rule = 2)
  # Modeling the variogram of the underlying Gaussian variable Y
  g.Yp <- vario.calc(db.data)
  g.Y <- g.Yp
  c.Yp <- var.Yp.mod - g.Yp@vardirs[[1]]@gg
  c.Y <- f.Y.vs.Yp(c.Yp)
  g.Y@vardirs[[1]]@gg <- 1 - c.Y
  g.Y@vars <- 1
  g.Y.mod <- model.auto(g.Y)
  # Define the grid on where simulation are performed
  db.grid <-
    db.grid.init(
      db.data,origin = c(min(pol$x),min(pol$y)),nodes = res,extend = c(max(pol$x) -
                                                                         min(pol$x),max(pol$y) - min(pol$y))
    )
  db.grid <- db.polygon(db.grid,db.poly)
  
  neigh <- neigh.read("neigh_pol.txt")
  # neigh <- neigh.input();
  # 2
  # 5
  # 100
  # n
  # n
  # 60
  
  # Define the interval limits for the gibbs sampling
  attach(db.data@items)
  Ymax <- Yp
  Ymin <- Yp
  Ymin[Ymin <= ycut] <- -10
  detach(db.data@items)
  db.data <-
    db.add(db.data,Ymax,type.locate = FALSE,loctype = "upper")
  db.data <-
    db.add(db.data,Ymin,type.locate = FALSE,loctype = "lower")
  db.data <- db.locate(db.data,"Ymax","upper")
  db.data <- db.locate(db.data,"Ymin","lower")
  rm(Ymin,Ymax)
  # Simulating gaussian values below ycut
  db.data <-
    gibbs(
      db = db.data,model = g.Y.mod,seed = as.numeric(Sys.time()),
      nbsimu = 1,rank.grf = 1,
      nboot = 100,niter = 1,percent = 0,toleps = 1,radix = NA,modify.target =
        T
    )
  db.data <- db.rename(db.data,"G1","Y")
  # conditional simulation of the gaussian variable
  db.simu <-
    simtub(
      dbin = db.data,dbout = db.grid,model = g.Y.mod,neigh = neigh,nbsimu = 1,seed = 
        as.numeric(Sys.time())
    )
  db.orig <- db.simu
  
  for (tspac in 1:length(transpacmulti)) {
    tmulti <- transpacmulti[tspac]
    cat(paste(
      "year: ",as.character(year),"    spacing multiplier: ",as.character(tmulti),"\n"
    ))
    #get sample estimates for given transect design
    for (t in 1:nbtransimu) {
      cat(paste("iteration:" ,as.character(t)),"\n")
      db.simu <- db.orig
      db.simu <- db.rename(db.simu,"Simu.Y.S1","Simu.Y")
      # transform the gaussian conditional simulation into raw conditional simulation
      db.simu <-
        db.add(db.simu,Simu.Z = f.phi(db.simu@items$Simu.Y,ycut,fct = f))
      # Visualization of the conditional simulation of Y
      db.simu <- db.locate(db.simu,"Simu.Y","z")
      db.simu <- db.locate(db.simu,seq(2,3),"x")
          plot(db.simu,col=tim.colors(1000),xlab="",ylab="",pos.legend=7,asp=1,xlim=c(min(pol$x)*1.2,max(pol$x)*2.5))
          plot(db.poly,add=TRUE,col=1,lwd=2)
          polygon(europa,col=8);box()
      # Visualization of the conditional simulation of Z > 0
      db.simu <- db.locate(db.simu,"Simu.Z","z")
      db.simu <- db.locate(db.simu,seq(2,3),"x")
          plot(db.simu,col=c(0,tim.colors(1000)),xlab="",ylab="",pos.legend=7,asp=1,xlim=c(min(pol$x)*1.2,max(pol$x)*2.5))
          plot(db.poly,add=TRUE,col=1,lwd=2)
          polygon(europa,col=8);box()
          points(txectdata$lon_dec,txectdata$lat_dec,type="l",col=2)
      
      simu.mean.sa <- mean(db.extract(db.simu,"Simu.Z"))
      
      #find points covered by transects and plot samples

### HORIZONTAL TXSECT OPTION ###
### HORIZONTAL TXSECT OPTION ###
          res <- seq(0,0.5,len=intperlong)[2]-seq(0,0.5,len=intperlong)[1]
          
          round.simu <-
            round(unique(db.extract(db.simu,"x2")) / res) * res
          
          ### randomise transect latitudes
          txectlats <- unique(txectdata$lat_dec)
          minpolylat <- min(pol$y) #minimum latitude of polygon
          maxpolylat <- max(pol$y) #maximum latitude of polygon
          txectspac <- (txectlats[2] - txectlats[1]) * tmulti
          txectstart <- runif(1,minpolylat,minpolylat + txectspac)
          randtxectlats <- seq(txectstart,maxpolylat,txectspac)
          ###
          
          round.txect <- round(randtxectlats / res) * res
          simu.latsamples <-
            unique(db.extract(db.simu,"x2"))[which(round.simu %in% round.txect)]
          
          simu.Zsamples <-
            db.extract(db.simu,"Simu.Z")[which(db.extract(db.simu,"x2") %in% simu.latsamples)]
          simu.x1samples <-
            db.extract(db.simu,"x1")[which(db.extract(db.simu,"x2") %in% simu.latsamples)]
          simu.x2samples <-
            db.extract(db.simu,"x2")[which(db.extract(db.simu,"x2") %in% simu.latsamples)]
          sample.data <-
            data.frame(long = simu.x1samples,lat = simu.x2samples,sa = simu.Zsamples)
          db.simusample <- db.create(sample.data)
          db.simusample <- db.locate(db.simusample,seq(2,3),"x")
          db.simusample <- db.locate(db.simusample,"sa","z")
                          plot(db.simusample,inches=3,asp=1/cos(mean(data$lat)*pi/180),title="Acoustic backscatter samples \nfrom conditonal simulations")
                          plot(db.poly,add=TRUE,col=1,lwd=2)
                          polygon(europa,col=8);box()
          
          sample.data <-
            transform(sample.data,txect.id = as.numeric(factor(sample.data$lat)))

### HORIZONTAL TXSECT OPTION ###
### HORIZONTAL TXSECT OPTION ###

# ### VERTICAL TXSECT OPTION ###
# ### VERTICAL TXSECT OPTION ###
#           res <- seq(0,1,len=intperlong)[2]-seq(0,1,len=intperlong)[1]
#           
#           round.simu <-
#             round(unique(db.extract(db.simu,"x1")) / res) * res
#           
#           ### randomise transect latitudes
#           txectlats <- unique(txectdata$lat_dec)
#           minpolylat <- min(pol$x) #minimum latitude of polygon
#           maxpolylat <- max(pol$x) #maximum latitude of polygon
#           txectspac <- (txectlats[2] - txectlats[1]) * tmulti * 2
#           txectstart <- runif(1,minpolylat,minpolylat + txectspac)
#           randtxectlats <- seq(txectstart,maxpolylat,txectspac)
#           ###
#           
#           round.txect <- round(randtxectlats / res) * res
#           simu.latsamples <-
#             unique(db.extract(db.simu,"x1"))[which(round.simu %in% round.txect)]
#           
#           simu.Zsamples <-
#             db.extract(db.simu,"Simu.Z")[which(db.extract(db.simu,"x1") %in% simu.latsamples)]
#           simu.x1samples <-
#             db.extract(db.simu,"x1")[which(db.extract(db.simu,"x1") %in% simu.latsamples)]
#           simu.x2samples <-
#             db.extract(db.simu,"x2")[which(db.extract(db.simu,"x1") %in% simu.latsamples)]
#           sample.data <-
#             data.frame(long = simu.x1samples,lat = simu.x2samples,sa = simu.Zsamples)
#           db.simusample <- db.create(sample.data)
#           db.simusample <- db.locate(db.simusample,seq(2,3),"x")
#           db.simusample <- db.locate(db.simusample,"sa","z")
#                           plot(db.simusample,inches=3,asp=1/cos(mean(data$lat)*pi/180),title="Acoustic backscatter samples \nfrom conditonal simulations")
#                           plot(db.poly,add=TRUE,col=1,lwd=2)
#                           polygon(europa,col=8);box()
#           
#           sample.data <-
#             transform(sample.data,txect.id = as.numeric(factor(sample.data$long)))          
# ### VERTICAL TXSECT OPTION ###
# ### VERTICAL TXSECT OPTION ###
          
      sample.mean.sa <-
        weighted.mean(
          aggregate(sample.data[, 3], list(sample.data$txect.id), mean)[,2],aggregate(sample.data[, 3], list(sample.data$txect.id), length)[,2]
        )
      sample.var.sa <-
        wtd.var(
          aggregate(sample.data[, 3], list(sample.data$txect.id), mean)[,2],aggregate(sample.data[, 3], list(sample.data$txect.id), length)[,2]
        )
      sample.se.sa <-
        sqrt(sample.var.sa) / sqrt(length(aggregate(
          sample.data[, 3], list(sample.data$txect.id), length
        )[,2]))
      sample.RSE.sa <- sample.se.sa / sample.mean.sa * 100
      
      yrRSE[t,1,tspac,y] <- simu.mean.sa
      yrRSE[t,2,tspac,y] <- sample.mean.sa
      yrRSE[t,3,tspac,y] <- sample.RSE.sa
      
    }
    
  }
  
  save(db.simu,file = file.path(uncertainty.dir,paste("simudata_stratum",as.character(stratum),"_year",as.character(year),".RData",sep = "")))
  
}

save(yrRSE,file = file.path(uncertainty.dir,paste("RSE_stratum",as.character(stratum),".RData",sep = "")))

end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste("COMPLETE IN",sprintf("%0.1f",round(time.taken)),"mins"))

}