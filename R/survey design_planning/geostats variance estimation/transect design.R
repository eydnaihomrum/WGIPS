#-------------------------------------------------------------------------------
# Code to format strata polygon files and transect design outputs using the
# survey design software 'Distance 6.2' (http://distancesampling.org/)
#
# download:
# http://distancesampling.org/Distance/old-versions/downloads/d62setup.exe
#
# By: Sascha F?ssler
# IMARES
# 9 April 2016
#-------------------------------------------------------------------------------

rm(list = ls());
library(rgdal)
library(raster)
path <-
  "~/Dropbox/survey demo/survey design"
try(setwd(path),silent = TRUE)
flatfile.dir <- file.path(".","Strata","2016","Flat File")
polygons.dir <- file.path(".","polygons")
transects.dir <- file.path(".","transects")
europa <- read.table("europa.txt",header = T)

# create polygon shape txt files with lon & lat boundaries from csv flat file
flatfile <-
  read.csv(file = file.path(flatfile.dir,"New HERAS Strata flat file.csv"),head =
             T)[-1]
strata <- unique(flatfile$ID)
stratasplit <- split(flatfile,flatfile$ID)
for (f in 1:length(strata)) {
  dat <-
    data.frame(lon = as.data.frame(stratasplit[f])[,1],lat = as.data.frame(stratasplit[f])[,2])
  write.table(
    dat,file = file.path(
      polygons.dir,paste("stratum_",strata[f],".LatLon.txt",sep = "")
    ),row.names = F,col.names = F,sep = "\t"
  )
  png(
    file.path(polygons.dir,paste("stratum ",strata[f],".png",sep = "")),units = "px", height =
      500,width = 500, bg = "white"
  )
  plot(
    dat$lon,dat$lat,type = "l",main = paste("stratum",strata[f]),lwd = 2,xlim =
      c(-4,13),ylim = c(51,62),col = 2,xlab = "longitude",ylab = "latitude"
  ); polygon(europa,col = 8);box()
  dev.off()
}

# load polygon shape txt files with lon & lat boundaries and convert to meter coordinates
strata.select <- as.character(c(41)) #select strata
polyfiles.latlon <-
  list.files(polygons.dir,pattern = paste(paste("_",strata.select,".LatLon",sep =
                                                  ""), collapse = "|"))
for (f in 1:length(polyfiles.latlon)) {
  polydata <-
    read.table(
      file = file.path(polygons.dir,polyfiles.latlon[f]),sep = "\t",head = F,col.names = c("lon","lat")
    )
  cord.dec <-
    SpatialPoints(cbind(polydata$lon, polydata$lat), proj4string = CRS("+proj=longlat"))
  # plot(data.frame(cord.dec)$coords.x1,data.frame(cord.dec)$coords.x2,type="l")
  cord.UTM <-
    spTransform(cord.dec, CRS("+init=epsg:4087")) #World Equidistant Cylindrical=4087
  d.converted <-
    data.frame(lon = cord.UTM$coords.x1/2,lat = cord.UTM$coords.x2)
  plot(d.converted$lon,d.converted$lat,type="l")
  write.table(
    d.converted,file = file.path(polygons.dir,paste(unlist(
      strsplit(polyfiles.latlon[f],"LatLon.txt")
    )[1],"Meters.txt",sep = "")),row.names = F,col.names = F,sep = "\t"
  )
}

#--- use these files in 'Distance 6.2' to create systematic random survey transects per stratum

# get transect information from 'Distance 6.2' survey design output (copied to txt file called 'stratum xx DistanceOutput.txt')
txectfiles.latlon <-
  list.files(transects.dir,pattern = "DistanceOutput.txt")
for (f in 1:length(txectfiles.latlon)) {
  txectdata <-
    read.table(file = file.path(transects.dir,txectfiles.latlon[f]),sep = "\t")
  txectdata <- txectdata[6:nrow(txectdata),] #exclude header
  transects <-
    as.numeric(gsub("[^0-9]","",as.character(txectdata[seq(1,length(txectdata),4)])))
  wp.data <-
    txectdata[sort(c(seq(2,length(txectdata),4),seq(3,length(txectdata),4)))]
  latlon <- as.numeric(unlist(strsplit(as.character(wp.data)," ")))
  latlon <- latlon[!is.na(latlon)]
  wp.lon <- latlon[seq(1,length(latlon),2)]*2
  wp.lat <- latlon[seq(2,length(latlon),2)]
  cord.UTM <-
    SpatialPoints(cbind(wp.lon, wp.lat), proj4string = CRS("+init=epsg:4087"))
  cord.dec <-
    spTransform(cord.UTM, CRS("+proj=longlat"))
  dec <- data.frame(coordinates(cord.dec))
  # plot(dec$wp.lon,dec$wp.lat,type = "l"); polygon(europa,col = 8);box()
  txect.df <-
    data.frame(
      lon.m = wp.lon,lat.m = wp.lat, lon.dec = dec$wp.lon, lat.dec = dec$wp.lat, transect = sort(rep(transects,2)), stratum = as.numeric(unlist(strsplit(
        txectfiles.latlon[f]," "
      ))[2])
    )
  evens <- txect.df[txect.df$transect %in% transects[transects %% 2 == 0],]
  odds <- txect.df[txect.df$transect %in% transects[transects %% 2 == 1],]
  turn.evens <- txect.df[rev(as.numeric(rownames(evens))),]
  turn.evens <- turn.evens[order(turn.evens$transect),]
  txect.df <- rbind(turn.evens,odds)
  txect.df <- txect.df[order(txect.df$transect),]
  plot(txect.df$lon.dec,txect.df$lat.dec,type="l"); polygon(europa,col = 8);box()
  write.table(
    txect.df,file = file.path(transects.dir,paste(unlist(
      strsplit(txectfiles.latlon[f],"DistanceOutput.txt")
    )[1],"waypoints.txt",sep = "")),row.names = F,col.names = F,sep = "\t"
  )
}
