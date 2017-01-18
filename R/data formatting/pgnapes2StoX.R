################################################################################
# PGNAPES (csv files) to StoX (xml files) conversion
#
# $Date: 2015-09-02 22:21:57 CEST $
#
# Author: Sascha Fässler
#
################################################################################
rm(list=ls())

pgnapes.biology <- read.csv("2014TRI2HERAS_Biology.csv",head=TRUE)
pgnapes.acoustic <- read.csv("2014TRI2HERAS_Acoustic.csv",head=TRUE)
pgnapes.acousticvalues <- read.csv("2014TRI2HERAS_AcousticValues.csv",head=TRUE)
pgnapes.logbook <- read.csv("2014TRI2HERAS_Logbook.csv",head=TRUE)
pgnapes.catch <- read.csv("2014TRI2HERAS_Catch.csv",head=TRUE)

### load required packages
require(XML)

### combine acoustic pgnapes data tables
acoustic.data <- merge(pgnapes.acoustic,pgnapes.acousticvalues,by=(names(pgnapes.acoustic)[1:7]),sort=TRUE)
acoustic.data[,15][is.na(acoustic.data[,15])]<-unique(pgnapes.acousticvalues[,8])
acoustic.data[,16][is.na(acoustic.data[,16])]<-min(pgnapes.acousticvalues[,9])
acoustic.data[,17][is.na(acoustic.data[,17])]<-max(pgnapes.acousticvalues[,10])
acoustic.data[,18][is.na(acoustic.data[,18])]<-0
acoustic.data <- setNames(acoustic.data,toupper(names(acoustic.data)))

### convert acoustic tables into StoX input format (XML)
len <- nrow(acoustic.data)
acoustic <- data.frame("nation"=numeric(len),"platform"=numeric(len),"cruise"=character(len),"start_time"=character(len),"log_start"=numeric(len),"integrator_dist"=numeric(len),"lat_start"=numeric(len),"lon_start"=numeric(len),"freq"=numeric(len),"transceiver"=numeric(len),"threshold"=numeric(len),"upper_interpret_depth"=numeric(len),"upper_integrator_depth"=numeric(len),"num_pel_ch"=numeric(len),"type"=character(len),"acocat"=numeric(len),"sa"=numeric(len),"ch"=numeric(len),"pel_ch_thickness"=numeric(len))

acoustic$nation <- 64
acoustic$platform <- 7156
acoustic$cruise <- acoustic.data$CRUISE
date <- format(as.Date(paste(acoustic.data$YEAR,acoustic.data$MONTH,acoustic.data$DAY,sep="-"), "%Y-%m-%d"),"%Y-%m-%d")
time <- format(as.POSIXlt(strptime(paste(acoustic.data$HOUR,acoustic.data$MIN,"00",sep=":"),format="%H:%M:%S")),"%H:%M:%S")
acoustic$start_time <- paste(date,time)
acoustic$log_start <- acoustic.data$LOG
acoustic$integrator_dist <- acoustic.data$LOGINT
acoustic$pel_ch_thickness <- acoustic.data$CHLOWDEPTH-acoustic.data$CHUPPDEPTH
acoustic$lat_start <- acoustic.data$ACLAT
acoustic$lon_start <- acoustic.data$ACLON
acoustic$freq <- acoustic.data$FREQUENCY * 1000
acoustic$transceiver <- 1	
acoustic$threshold <- acoustic.data$SV_THRESHOLD	
acoustic$num_pel_ch <- 1 #workaround for now, but needs to be max channel number per interval: max(acoustic.data$CHLOWDEPTH)  
acoustic$upper_interpret_depth <- 0	
acoustic$upper_integrator_depth <- 0	
acoustic$type  <- "P"
acoustic$acocat[which(acoustic.data$SPECIES=="WHB")] <- 24 # BLUE WHITING
acoustic$acocat[which(acoustic.data$SPECIES=="HER")] <- 12 # HERRING
acoustic$acocat[which(acoustic.data$SPECIES=="SPR")] <- 13 # SPRAT???
acoustic$sa <- acoustic.data$SA	
acoustic$ch <- 1 #only one channel for now... otherwise: acoustic.data$CHLOWDEPTH/unique(acoustic.data$CHLOWDEPTH-acoustic.data$CHUPPDEPTH)	

### order by log
acoustic <- acoustic[order(acoustic$log_start),]

### convert to XML file
acoustic.xml = newXMLNode("echosounder_dataset")
newXMLNode(names(acoustic)[1], unique(acoustic[1]), parent = acoustic.xml)
newXMLNode(names(acoustic)[2], unique(acoustic[2]), parent = acoustic.xml)
newXMLNode(names(acoustic)[3], unique(acoustic[3]), parent = acoustic.xml)
for (i in 1:length(unique(acoustic$log_start))){
o = newXMLNode("distance_list", parent = acoustic.xml)
  p = newXMLNode("distance", attrs = c(as.vector((acoustic)[i,4:5])), parent = acoustic.xml)
    q = newXMLNode("frequency", attrs = c(as.vector((acoustic)[i,9:10])))
      r = newXMLNode("ch_type", attrs = c(type = acoustic[i,15]), parent = acoustic.xml)
        s = newXMLNode("sa_by_acocat", attrs = c(acocat = acoustic[i,16]), parent = acoustic.xml)
          sa = newXMLNode(names(acoustic)[17],acoustic[i,17])
          addAttributes(sa, ch=acoustic[i,18])
          addChildren(s,sa)
      addChildren(r,s)
    frequencies = c(newXMLNode(names(acoustic)[11],acoustic[i,11]),newXMLNode(names(acoustic)[14],acoustic[i,14]),newXMLNode(names(acoustic)[12],acoustic[i,12]),newXMLNode(names(acoustic)[13],acoustic[i,13]),r)
    addChildren(q, frequencies)  
  distances = c(newXMLNode(names(acoustic)[6], acoustic[i,6]),newXMLNode(names(acoustic)[19],acoustic[i,19]),newXMLNode(names(acoustic)[7],acoustic[i,7]),newXMLNode(names(acoustic)[8],acoustic[i,8]),q)
  addChildren(p, distances)
addChildren(o, p)
}

### save acoustic XML file
saveXML(acoustic.xml,file=paste("Echosounder-HERAS.xml",sep=""),prefix='<?xml version="1.0" encoding="UTF-8"?>')


### get required catch pgnapes data tables
catch.data <- pgnapes.catch[which(pgnapes.catch$Species=="HER"),]
logbook.data <- pgnapes.logbook[which(pgnapes.logbook$StType=="PTRAWL"),]
biology.data <- pgnapes.biology[which(pgnapes.biology$Species=="HER"),]

spptrawls <- unique(biology.data[,c(4,7)])$Station
logbook.data <- logbook.data[which(logbook.data$Station %in% spptrawls),]
catch.data <- catch.data[which(catch.data$Station %in% spptrawls),]

### convert biotic tables into StoX input format (XLM) [only use trawls containing target species]
len <- (length(spptrawls) * (20 + 6)) + #trawl header/closer 
        (length(spptrawls) * 6) + #spp samples header
          nrow(biology.data) #nr of samples
biotic <- data.frame("year"=numeric(len),"missiontype"=numeric(len),"missiontypename"=character(len),"callsignal"=character(len),"platform"=numeric(len),"missionnumber"=numeric(len),"cruise"=character(len),"serialno"=numeric(len),"nation"=numeric(len),"platform2"=numeric(len),"station"=numeric(len),"startdate"=character(len),"starttime"=character(len),"stoptime"=character(len),"latitudestart"=numeric(len),"longitudestart"=numeric(len),"system"=numeric(len),"area"=numeric(len),"location"=numeric(len),"bottomdepthstart"=numeric(len),"bottomdepthstop"=numeric(len),"fishingdepthmax"=numeric(len),"fishingdepthmin"=numeric(len),"gear"=numeric(len),"gearspeed"=numeric(len),"startlog"=numeric(len),"gearcondition"=numeric(len),"trawlquality"=numeric(len),"samplenumber"=numeric(len),"species"=numeric(len),"noname"=character(len),"aphia"=numeric(len),"sampletype"=numeric(len),"conservation"=numeric(len),"producttype"=numeric(len),"weight"=numeric(len),"lengthmeasurement"=character(len),"lengthsamplecount"=numeric(len),"specimenno"=numeric(len),"lengthunit"=numeric(len),"length"=numeric(len),"producttype2"=numeric(len),"weight3"=numeric(len),"sex"=numeric(len),"specialstage"=numeric(len),"stomachfillfield"=numeric(len),"no"=numeric(len),"age"=numeric(len),"readability"=numeric(len),"count"=numeric(len),"sampleproducttype"=numeric(len),"lengthsampleweight"=numeric(len),"specimensamplecount"=numeric(len),"agingstructure"=numeric(len),"wirelength"=numeric(len))

biotic$year <- as.numeric(unique(logbook.data$Year))
biotic$missiontype <- 4
biotic$missiontypename <- "research vessel"
biotic$callsignal <- as.character(unique(logbook.data$Vessel))
biotic$platform <- 7156
biotic$missionnumber <- as.numeric(paste(unique(logbook.data$Year),0,min(unique(logbook.data$Month)),sep=""))
biotic$cruise <- as.character(unique(logbook.data$Cruise))
biotic$serialno <- rep((1:length(unique(biology.data[,c(4,7)])$Station)),data.frame(table(biology.data$Station))$Freq + 32)
date <- format(as.Date(paste(logbook.data$Day,logbook.data$Month,logbook.data$Year,sep="/"), "%d/%m/%Y"),"%d/%m/%Y")
time <- format(as.POSIXlt(strptime(paste(logbook.data$Hour,logbook.data$Min,"00",sep=":"),format="%H:%M:%S")),"%H:%M:%S")
biotic$startdate <- as.character(biotic$startdate)
biotic$starttime <- as.character(biotic$starttime)
biotic$stoptime <- as.character(biotic$stoptime)
biotic$noname <- as.character(biotic$noname)
biotic$lengthmeasurement <- as.character(biotic$lengthmeasurement)
s <- 0
for(r in 1:length(unique(biology.data[,c(4,7)])$Station)){
  spprows <- data.frame(table(biology.data$Station))$Freq[r]
  trawlrows <- 32 + spprows
  biotic$nation[1+s] <- 64
  biotic$platform2[2+s] <- 7156
  biotic$station[3+s] <- unique(logbook.data$Station)[r]
  biotic$startdate[4+s] <- date[r]
  biotic$starttime[5+s] <- time[r]
  biotic$stoptime[6+s] <- time[r]
  biotic$latitudestart[7+s] <- logbook.data$Lat[r]
  biotic$longitudestart[8+s] <- logbook.data$Lon[r]
  biotic$system[9+s] <- 2
  biotic$area[10+s] <- 0
  biotic$location[11+s] <- 0
  biotic$bottomdepthstart[12+s] <- 0
  biotic$bottomdepthstop[13+s] <- 0
  biotic$fishingdepthmax[14+s] <- 0
  biotic$fishingdepthmin[15+s] <- 0
  biotic$gear[16+s] <- 3500
  biotic$gearspeed[17+s] <- 3.5
  biotic$startlog[18+s] <- logbook.data$Log[r]
  biotic$gearcondition[19+s] <- 1
  biotic$trawlquality[20+s] <- 2
  biotic[(21+s):(21+s+spprows+6),29] <- 1
  biotic[(21+s):(21+s+spprows+6),30] <- 161789 #sprat; 161722 #herring ;
  biotic[(21+s):(21+s+spprows+6),31] <- "Sprattus sprattus" #sprat; "Clupea harengus" #herring ;
  biotic[(21+s):(21+s+spprows+6),32] <- 126425 #sprat; 126417 #herring; 
  biotic$sampletype[21+s] <- 20
  biotic$conservation[22+s] <- 1
  biotic$producttype[23+s] <- 1
  biotic$weight[24+s] <- catch.data$Catch[r]
  biotic$lengthmeasurement[25+s] <- "E"
  biotic$lengthsamplecount[26+s] <- spprows
  biotic[(27+s):(27+s+spprows-1),39] <- biology.data$Recnr[biology.data$Station == unique(logbook.data$Station)[r]]
  biotic[(27+s):(27+s+spprows-1),40] <- 2
  biotic[(27+s):(27+s+spprows-1),41] <- biology.data$Length[biology.data$Station == unique(logbook.data$Station)[r]]/1e2
  biotic[(27+s):(27+s+spprows-1),43] <- biology.data$Weight[biology.data$Station == unique(logbook.data$Station)[r]]/1e3
  biotic[(27+s):(27+s+spprows-1),44] <- abs(biology.data$Sex[biology.data$Station == unique(logbook.data$Station)[r]]-3) # 1=female, 2=male... not other way round
  biotic[(27+s):(27+s+spprows-1),45] <- biology.data$Maturation[biology.data$Station == unique(logbook.data$Station)[r]] # 2=immature, 4=matture in the NL case...
  biotic[(27+s):(27+s+spprows-1),48] <- biology.data$AgeOtolith[biology.data$Station == unique(logbook.data$Station)[r]]
  biotic$producttype2[which(biotic$weight3 > 0)] <- 1
  biotic$no[which(biotic$weight3 > 0)] <- 1
  biotic$readability[which(biotic$weight3 > 0)] <- 2
  biotic$count[27+s+spprows] <- 1000 #constant for now... raised length sample count
  biotic$sampleproducttype[28+s+spprows] <- 1
  biotic$lengthsampleweight[29+s+spprows] <- 10 #constant for now... raised length sample weight
  biotic$specimensamplecount[30+s+spprows] <- length(na.omit(biology.data$Weight[biology.data$Station == unique(logbook.data$Station)[r]]))
  biotic$agingstructure[31+s+spprows] <- 2
  biotic$wirelength[32+s+spprows] <- 350 # constant for now
  s <- s + trawlrows
}

### convert to XML file
biotic.xml = newXMLNode("mission",attrs = c(as.vector(unique(biotic[1:6]))))
newXMLNode(names(biotic)[7], unique(biotic[7]),parent=biotic.xml) #cruise
for (i in 1:length(unique(logbook.data$Station))){
  trdat <- biotic[which(biotic$serialno == i),]
  lengthsamples <- trdat[which(trdat$length>0 & is.na(trdat$weight3)),]
  weightsamples <- trdat[which(trdat$weight3>0),]
  o = newXMLNode("fishstation", attrs = c(unique(as.vector((trdat)[8]))),parent = biotic.xml)
  p = newXMLNode("catchsample",attrs = c(as.vector(unique(trdat[which(trdat[29]>0),29:32]))))
  catchsamples = c(newXMLNode(names(trdat)[33],max(trdat[,33])),
                   newXMLNode(names(trdat)[34],max(trdat[,34])),
                   newXMLNode(names(trdat)[35],max(trdat[,35])),
                   newXMLNode(names(trdat)[36],max(trdat[,36])),
                   newXMLNode(names(trdat)[50],max(trdat[,50])),
                   newXMLNode(names(trdat)[51],max(trdat[,51])),
                   newXMLNode(names(trdat)[37],max(trdat[,37])),
                   newXMLNode(names(trdat)[52],max(trdat[,52])),
                   newXMLNode(names(trdat)[38],max(trdat[,38])),
                   newXMLNode(names(trdat)[53],max(trdat[,53])),
                   newXMLNode(names(trdat)[54],max(trdat[,54])))
  addChildren(p, catchsamples)  
  for(ls in 1:nrow(lengthsamples)){
    q = newXMLNode("individual", attrs = c(specimenno = ls),parent=catchsamples)
    individuals = c(newXMLNode(names(lengthsamples)[40],lengthsamples[ls,40]),
                    newXMLNode(names(lengthsamples)[41],lengthsamples[ls,41]))
    addChildren(q,individuals)
    addChildren(p,q)
  }
  for(ws in 1:nrow(weightsamples)){
    q = newXMLNode("individual", attrs = c(specimenno = nrow(lengthsamples)+ws),parent=catchsamples)
    r = newXMLNode("agedetermination", attrs = c(unique(weightsamples[47])))
    age = c(newXMLNode(names(weightsamples)[48],weightsamples[ws,48]))
    individuals = c(newXMLNode("producttype",weightsamples[ws,42]),
                    newXMLNode("weight",weightsamples[ws,43]),
                    newXMLNode(names(weightsamples)[40],weightsamples[ws,40]),
                    newXMLNode(names(weightsamples)[41],weightsamples[ws,41]),
                    newXMLNode(names(weightsamples)[44],weightsamples[ws,44]),
                    newXMLNode(names(weightsamples)[45],weightsamples[ws,45]),
                    newXMLNode(names(weightsamples)[46],weightsamples[ws,46]),r)
    addChildren(r,age)
    addChildren(q,individuals)
    addChildren(p,q)
  }
  fishstation = c(newXMLNode(names(trdat)[9],max(trdat[,9])),
                  newXMLNode(names(trdat)[10],max(trdat[,10])),
                  newXMLNode(names(trdat)[11],max(trdat[,11])),
                  newXMLNode(names(trdat)[12],max(trdat[,12])),
                  newXMLNode(names(trdat)[13],max(trdat[,13])),
                  newXMLNode(names(trdat)[14],max(trdat[,14])),
                  newXMLNode(names(trdat)[15],max(trdat[,15])),
                  newXMLNode(names(trdat)[16],sum(trdat[,16])),
                  newXMLNode(names(trdat)[17],max(trdat[,17])),
                  newXMLNode(names(trdat)[18],max(trdat[,18])),
                  newXMLNode(names(trdat)[19],max(trdat[,19])),
                  newXMLNode(names(trdat)[20],max(trdat[,20])),
                  newXMLNode(names(trdat)[21],max(trdat[,21])),
                  newXMLNode(names(trdat)[22],max(trdat[,22])),
                  newXMLNode(names(trdat)[23],max(trdat[,23])),
                  newXMLNode(names(trdat)[24],max(trdat[,24])),
                  newXMLNode(names(trdat)[25],max(trdat[,25])),
                  newXMLNode(names(trdat)[26],max(trdat[,26])),
                  newXMLNode(names(trdat)[27],max(trdat[,27])),
                  newXMLNode(names(trdat)[28],max(trdat[,28])),p)
  addChildren(o,fishstation)  
}

### save biology XML file
saveXML(biotic.xml,file=paste("Biotic-HERAS.xml",sep=""),prefix='<?xml version="1.0" encoding="UTF-8"?>')
