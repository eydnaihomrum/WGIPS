paralltrx <- function(inndata, rev.entrance=F) {

# paralltrx
#paralltr <- function(inndata){

datapath <- "D:/Sascha/Projects/WGIPS/WGIPS 2016/HERAS post-cruise/Harbitz and Johnsen II/Kartdata/Toktomrader2015/"
fieldname <- read.table(paste(datapath,"fieldname.txt", sep="") ,sep = '\n')
fieldnameText <- read.table(paste(datapath,"fieldnamesText.txt", sep=""))
lonlat <- read.table(paste(datapath,"lonlatx.txt", sep=""))
lonlat12 <- read.table(paste(datapath,"lonlat12x.txt", sep=""))

library("splancs")
library("sp")
library(geosphere)

# INTERACTIVE INPUT
print(fieldname)
no <- readline('Field number             no = ? ')
t <- readline('time available in hours,  t = ? ')
vel <- readline('vessel speed in knots,  vel = ? ')
nsim <- readline('no of runs      ,      nsim = ? ')
no <- as.numeric(no)
t <- as.numeric(t)
vel <- as.numeric(vel)
nsim <- as.numeric(nsim)
name.field <- as.character(fieldnameText[fieldnameText[,1]==no,2])
print(name.field)
if (nsim == 1) {
	nsim == 0
	}

# EXTRACTION AND FORMATTING OF POLYGON NODES AND ENTRANCE AND EXIT POINTS
lbno <- lonlat[which(lonlat[,3] == no),1:2]			# picks polygon nodes for field no
lon <- lbno[,1]											# polygon node decimal longitudes
lat <- lbno[,2]											# polygon node decimal latitudes
lb12 <- lonlat12[which(lonlat12[,5] == no),1:4]	# picks entrance and exit points for field no
if(rev.entrance==T) lb12 <- lonlat12[which(lonlat12[,5] == no),c(3,4,1,2)] # reverse entrance and exit
lb12 <- as.numeric(lb12)								# converts to numeric (lb12 is one line)
lon1 <- lb12[1]				# decimal longitude of entrance point
lat1 <- lb12[2]				# decimal latitude of entrance point
lon2 <- lb12[3]				# decimal longitude of exit point
lat2 <- lb12[4]				# decimal latitude of exit point

n <- length(lon);			# number of nodes in polygon 
lon <- matrix(lon,n,1)		# ascertains lon to be column vector
lat <- matrix(lat,n,1)		# ascertains lat to be column vector
# ascertains lon(1) = lon(end) and lat(1) = lat(end) if the difference of entrance and exit points
# exceed 1e-6 deg for either lon-values, lat-values or both
if (abs(lon[1]-lon[n]) > 1e-6 | abs(lat[1]-lat[n]) > 1e-6)	# TRUE if difference exceeds threshold
	{
    lon <- c(lon, lon[1]);	# adds a new last lon-value identical to first one
    lat <- c(lat, lat[1]);	# adds a new last lat-value identical to first one
    n <- n + 1;				# adjusts n to define the length of lon and lat
    }
lon <- matrix(lon,n,1)		# ascertains lon to be a vertical vector
lat <- matrix(lat,n,1)		# ascertains lat to be a vertical vector

inndata <- cbind(lon,lat,lon1,lat1,lon2,lat2)
#	browser()
	y <- subconvex(inndata)
	uc <- y$uc
	vc <- y$vc
#	vel <- inndata[1,7]
#	t <- inndata[1,8]
	c <- ccalc(uc,vc,vel,t)
#	fac <- inndata[1,9]
#	if (fac == 0) {
	if (nsim == 1) {
		plot(lon,lat,type='l')
		}
	for (jsim in 1:max(1,nsim)) {
		fac <- runif(1)
#		}

	ur1 <- fac*c
	urnd <- seq(min(uc)+ur1, max(uc),by = c)
	nurnd <- length(urnd)
	uupp <- y$uupp
	ulow <- y$ulow
	vupp <- y$vupp
	vlow <- y$vlow
	q <- subyz(urnd, nurnd, uupp, ulow, vupp, vlow)
	yzu <- q$yz[,1]
	yzd <- q$yz[,2]
	A <- y$A
	c0 <- y$c0
	lon0 <- y$lon0
	lat0 <- y$lat0
	zig <- subuplow(urnd,nurnd,uupp,vupp,ulow,vlow,yzu,yzd,A,c0,lon0,lat0)
	lines(zig$ziglon1,zig$ziglat1,col = 'red')
	}
	zig <- subuplow(urnd,nurnd,uupp,vupp,ulow,vlow,yzu,yzd,A,c0,lon0,lat0)
                       
	lonrnd <- zig$ziglon1
	latrnd <- zig$ziglat1
  
  dist1 <- 0
  for(i in 1:(length(lonrnd)-1)){
    dist1 <- dist1 + (spDistsN1(cbind(lonrnd,latrnd),cbind(lonrnd,latrnd)[i,], longlat=T)[i+1])
    }
  dist1 <- dist1/1.852  

  #dist1 <- spDistsN1(cbind(lonrnd,latrnd),cbind(lonrnd,latrnd)[1,], longlat=T)[i+1]
  out <- list(name.field=name.field,
        cruiseline = as.data.frame(cbind(lonrnd,latrnd)), 
#        Area = A,
        Distance = dist1
#        Sur.cov = dist1/ sqrt(A)
        ) 

  ## Create file to MaxSea
  ## Create ascii fil til import i MaxSea (se manual for code)
  item.type <- 257  # Linje
  item.id <- 20      # Setter id-coder på linjene
  item.col <- 1     # Blå
  north <- "N"
  west <- "E"


  ##
  filnavn.ut <- paste("D:/Sascha/Projects/WGIPS/WGIPS 2016/HERAS post-cruise/Harbitz and Johnsen II/Kartdata/rutelinjer/",name.field,".asc",sep="")
  ant.dec.pos <- 3
  n.row <- nrow(out$cruiseline)
                      
  out1 <- cbind(rep(item.type,n.row),
            rep(item.id,n.row),
            rep(item.col,n.row),
            round(out$cruiseline$lat,ant.dec.pos),
            rep(north,n.row),
            round(out$cruiseline$lonrnd,ant.dec.pos),
            rep(west,n.row))
            
  write.table(out1,file=filnavn.ut,row.names=F, col.names=F,sep=",",quote=F)

  polyg <- cbind(lon,lat)
  polyg.sp = SpatialPolygons(list(Polygons(list(Polygon(polyg)), "x")))
  #plot(polyg.sp, main=name.field)
  #plot(lon,lat,type = 'l',main=name.field)
  points(c(lon1,lon2),c(lat1,lat2),col=2)
  text(c(lon1,lon2),c(lat1,lat2),labels=1:2,col=2)
  lines(lonrnd, zig$ziglat1, col="red")
  #points(spsample(polyg.sp, n = 10, "random"), pch = 2, col=2)
  #spsample(polyg.sp, n = 10, "random")

	#out <- list(lon=lon,lat=lat,zig)
	#track.description <- paste("D:/Sascha/Projects/WGIPS/WGIPS 2016/HERAS post-cruise/Harbitz and Johnsen II/Kartdata/rutelinjer/",name.field,".asc",sep="")
	
  tmp <- as.data.frame(cbind(floor(out$cr$lonrnd), round(60*(out$cr$lonrnd- floor(out$cr$lonrnd)),3), 
    floor(out$cr$latrnd), round(60*(out$cr$latrnd- floor(out$cr$latrnd)),3)))
  names(tmp) <- c("LonDeg","LonMin","LatDeg","LatMin")  
  name.field
  vel # Speed
  dist1 ## Sailing distance
  area.nm2 <- areaPolygon(polyg)/1852/1852 # Meters2
  Sur.cov <- dist1/ sqrt(area.nm2)
  print(tmp)  
  ## TABLE
  filnavn1 <- paste("D:/Sascha/Projects/WGIPS/WGIPS 2016/HERAS post-cruise/Harbitz and Johnsen II/Kartdata/rutelinjer/INFO",name.field,".txt",sep="")
  capture.output( cat("\n",format(c(date()),width=20, justify = "left")), file=filnavn1)
  capture.output( cat("\n", "Stratum (Toktområde)    ", format(name.field,width=7, justify ="right")), file=filnavn1, append=T)
  capture.output( cat("\n", "Speed and time available", format(c(paste(vel,"knots"),paste(t, "h")),width=7, justify ="right")), file=filnavn1, append=T)
  capture.output( cat("\n", "Stratum area (n.mi2)    ", format(area.nm2,width=7, justify ="right")), file=filnavn1, append=T)
  capture.output( cat("\n", "Sailing distance (n.mi) ", format(dist1,width=7, justify ="right")), file=filnavn1, append=T)
  capture.output( cat("\n", "Survey coverage         ", format(Sur.cov,width=7, justify ="right")), file=filnavn1, append=T)
  capture.output( cat("\n", " "), file=filnavn1, append=T)
  capture.output( cat("\n", "Transect positions      "), file=filnavn1, append=T)
  capture.output( cat("\n", " "), file=filnavn1, append=T)
  capture.output(cat("\n", " "), tmp, file= filnavn1, append=T)  
 # browser()
  filnavn2 <- paste("D:/Sascha/Projects/WGIPS/WGIPS 2016/HERAS post-cruise/Harbitz and Johnsen II/Kartdata/rutelinjer/TRACK",name.field,".txt",sep="")
  out2 <- as.data.frame(cbind(1,out$cr))
  names(out2) <- c("Line","Longitude", "Latitude")
  write.csv(out2, file=filnavn2,row.names=F)
}


subconvex <- function(x) {
	lon <- x[,1]
	lat <- x[,2]
	lon1 <- x[1,3]
	lat1 <- x[1,4]
	lon2 <- x[1,5]
	lat2 <- x[1,6]
	lonc <- lon
	latc <- lat
	n = length(lonc)
  #browser()
	if (lonc[1] == lonc[n] & latc[1] == latc[n]) {
		lonc = lonc[2:n];
		latc = latc[2:n]
		}
	j <- which(diff(lonc) == 0 & diff(latc) == 0)
	nj <- length(j)
	if (nj > 0){
		ind = setdiff(1:length(lonc), j+1);
		lonc = lonc(ind);
		latc = latc(ind)
		}
	lon0 <- mean(lonc)
	lat0 <- mean(latc)
	c0 <- cos(lat0/180*pi);
	x1 <- (lon1-lon0)*60*c0
	y1 <- (lat1-lat0)*60
	x2 <- (lon2-lon0)*60*c0
	y2 <- (lat2-lat0)*60
	tet <- atan2(y2-y1,x2-x1)
	A = matrix(cbind(cos(tet), sin(tet), -sin(tet), cos(tet)),2,2,byrow = T)
	#out <- cbind(x1,y1)
	xc <- (lonc-lon0)*60*c0
	yc <- (latc-lat0)*60;
	uv <- A%*%t(cbind(xc,yc))
	uc <- uv[1,]
	vc <- uv[2,]
	n <- length(uc)
	jmin <- which(uc == min(uc))
	jmax <- which(uc == max(uc))
	if (length(jmin) == 2){
		j1 <- min(jmin);
		j2 <- max(jmin)
		if (j2 - j1 == 1) {
			uc <- t(c(uc[1:j1], uc[j1], uc[j2:n]));
			uc <- t(uc);
			vc <- t(c(vc[1:j1], (vc[j1]+vc[j2])/2, vc[j2:n]));
			vc <- t(vc);
			jmin <- j1[1] + 1
		}
		if (j2 - j1 > 1) {
			uc <- t(c(uc[j1], uc));
			uc <- t(uc);
			vc <- t(c((vc[j1]+vc[j2])/2, vc));
			vc <- t(vc);
			jmin <- 1
		}
	}
	if (length(jmax) == 2) {
		n <- length(uc);
		jmax <- which(uc == max(uc));
		j1 <- min(jmax);
		j2 <- max(jmax);
		if (j2 - j1 == 1){
			uc <- t(c(uc[1:j1], uc[j1], uc[j2:n]));
			uc <- t(uc);
			vc <- t(c(vc[1:j1], (vc[j1]+vc[j2])/2, vc[j2:n]));
			vc <- t(vc);
			jmax <- j1 + 1
		}
		if (j2 - j1 > 1) {
			uc <- t(c(uc[j1], uc));
			uc <- t(uc);
			vc <- t(c((vc[j1]+vc[j2])/2, vc));
			vc <- t(vc);
			jmax <- 1
		}
	}
	n <- length(uc)
	if (jmax > jmin){
		ind1 <- seq(jmin,jmax,1);
		ind2 <- c(seq(jmax,n,1), 1:jmin);
		nind2 <- length(ind2);
		ind2 <- ind2[nind2:1];
#		ind2 <- ind2[order(ind2, decreasing = T)]
	}
	if (jmax < jmin) {
		ind1 <- c(seq(jmin,n,1), 1:jmax);
#		ind2 <- jmax:jmin;
		ind2 <- jmin:jmax;
#		ind2 <- ind2[order(ind2,decreasing = T)]
	}
	if (mean(vc[ind1]) < mean(vc[ind2])){
		indupp <- ind2;
		indlow <- ind1
	}
	if (mean(vc[ind1]) > mean(vc[ind2])){
		indupp <- ind1;
		indlow <- ind2
	}
	dupp <- diff(uc[indupp])
	jupp <- which(dupp < 0)
	nupp <- length(jupp)
	if (nupp > 0){
		for (j in 1:nupp){
			j1 <- indupp[jupp[j]]
			j2 <- indupp[jupp[j]+1]
			if (vc[j1] < vc[j2]) {
				uc[j1] <- uc[j2]
			}
			if (vc[j1] > vc[j2]) {
				uc[j2] <- uc[j1]
			}
		}
	}
	dlow <- diff(uc[indlow])
	jlow <- which(dlow < 0)
	nlow <- length(jlow)
	if (nlow > 0){
		for (j in 1:nlow){
			j1 <- indlow[jlow[j]]
			j2 <- indlow[jlow[j] + 1]
			if (vc[j1] < vc[j2]) {
				uc[j1] <- uc[j2]
			}
			if (vc[j1] > vc[j2]) {
				uc[j2] <- uc[j1]
			}
		}
	}
	uupp <- uc[indupp]
	vupp <- vc[indupp]
	ulow <- uc[indlow]
	vlow <- vc[indlow]
	lonupp <- lonc[indupp]
	latupp <- latc[indupp]
	lonlow <- lonc[indlow]
	latlow <- latc[indlow]
	vm <- mean(vupp) - mean(vlow)
	out <- list(uc = uc, 
	vc = vc, 
	uupp = uupp,
	vupp = vupp,
	ulow = ulow,
	vlow =vlow,
	lonc = lonc,
	latc = latc,
	lon0 = lon0,
	lat0 = lat0,
	c0 = c0,
	vm = vm,
	A = A)
}

ccalc <- function(uc,vc,vel,t) {
	poly <- cbind(uc,vc)
	area <- areapl(poly)
	umin <- min(uc)
	umax <- max(uc)
	L <- vel*t
	du <- umax - umin
	c <- area/(vel*t - du)
	if (c < 0) {
		c <- 0.5*du
		}
		out <- c
	}
	
subyz <- function(urnd, nurnd, uupp, ulow, vupp, vlow) {
	yzu <- matrix(0,nurnd,1)
	yzd <- yzu
	for (j in 1:nurnd) {  # loop to identify (u,v) positions of nodes
    ju <- which(urnd[j]-uupp < 0);  # to identify actual upper polygon segment
    ju1 <- ju[1]-1;  # index to first point of upper polygon segment with node j
    ju2 <- ju1 + 1;  # index to second point of upperpolygon segment with node j
    yzu[j] <- vupp[ju1] + (urnd[j]-uupp[ju1])/(uupp[ju2]-uupp[ju1])*(vupp[ju2]-vupp[ju1]);
    # yzu = y-values of nodes in upper polygon
    jd <- which(urnd[j]-ulow < 0); # to identify actual lower polygon segment
    jd1 <- jd[1]-1;  # index to first point of lower polygon segment with node j
    jd2 <- jd1+1;  # index to second point of lower polygon segment with node j
    yzd[j] <- vlow[jd1] + (urnd[j]-ulow[jd1])/(ulow[jd2]-ulow[jd1])*(vlow[jd2]-vlow[jd1]);
    }
    # yzd = y-values of nodes in lower polygondasfasdf
    yz <- cbind(yzu,yzd)
    out <- list(yz = yz)
}
    
subuplow <- function(urnd,nurnd,uupp,vupp,ulow,vlow,yzu,yzd,A,c0,lon0,lat0) {
	#zigu1 <- NULL
	#zigv1 <- NULL
	#zigu2 <- NULL
	#zigv2 <- NULL
	for (j in 1:nurnd-1) {
		jodd <- j%%2
		ku <- which(uupp > urnd[j] & uupp  < urnd[j + 1])
		kd <- which(ulow > urnd[j] & ulow < urnd[j + 1])
		if (length(ku) > 0) {
			du1 <- c(urnd[j], uupp[ku], urnd[j+1]);
			dv1 <- c(yzu[j], vupp[ku], yzu[j+1])
			}
		if (length(ku) == 0) {
			du1 <- c(urnd[j], urnd[j+1]);
			dv1 <- c(yzu[j], yzu[j+1]	)
			}
		if (length(kd) > 0) {
			du2 <- c(urnd[j], ulow[kd], urnd[j+1]);
			dv2 <- c(yzd[j], vlow[kd], yzd[j+1])
			}
		if (length(kd) == 0) {
			du2 <- c(urnd[j], urnd[j+1]);
			dv2 <- c(yzd[j], yzd[j+1])
			}
		if (jodd == 1 & j == 1) {
			zigu1 <- du1;
			zigv1 <- dv1;
			zigu2 <- du2;
			zigv2 <- dv2
			}
		if (jodd == 0 & j == 1) {
			zigu1 <- du2;
			zigv1 <- dv2;
			zigu2 <- du1;
			zigv2 <- dv1
			}		
		if (jodd == 1 & j > 1) {
			zigu1 <- c(zigu1, du1);
			zigv1 <- c(zigv1, dv1);
			zigu2 <- c(zigu2, du2);
			zigv2 <- c(zigv2, dv2);
			}
    	if (jodd == 0 & j > 1) {
			zigu1 <- c(zigu1, du2); # updates successive u-values in transect one
        	zigv1 <- c(zigv1, dv2); # updates successive v-values in transect one
        	zigu2 <- c(zigu2, du1); # updates successive u-values in transect two
        	zigv2 <- c(zigv2, dv1); # updates successive v-values in transect two
        	}
		}
		zigu1 <- c(zigu2[1], zigu1)
		zigv1 <- c(zigv2[1], zigv1)
		zigu2 <- c(zigu1[2], zigu2)
		zigv2 <- c(zigv1[2], zigv2)
		
		zigu1 <- c(zigu1, zigu2[length(zigu2)])
		zigv1 <- c(zigv1, zigv2[length(zigv2)])
		zigu2 <- c(zigu2, zigu1[length(zigu1) - 1])
		zigv2 <- c(zigv2, zigv1[length(zigv1) - 1])       	
		
		j <- which(uupp < urnd[1])
		if (length(j) > 0) {
			zigu2 <- c(uupp[j], zigu2);
			zigv2 <- c(vupp[j], zigv2)
			}
		j <- which(ulow < urnd[1])
		if (length(j) > 0) {
			zigu1 <- c(ulow[j], zigu1);
			zigv1 <- c(vlow[j], zigv1)
			}
		ju <- which(uupp > urnd[length(urnd)])
		jd <- which(ulow > urnd[length(urnd)])
		jrnd <- nurnd%%2
		if (jrnd == 0) {
			if (length(ju) > 0) {
				zigu2 <- c(zigu2, uupp[ju]);
				zigv2 <- c(zigv2, vupp[ju])
				}
			if (length(jd) > 0) {
				zigu1 <- c(zigu1, ulow[jd])
				zigv1 <- c(zigv1, vlow[jd])
				}
			}
		if (jrnd == 1) {
			if (length(ju) > 0) {
				zigu1 <- c(zigu1, uupp[ju])
				zigv1 <- c(zigv1, vupp[ju])
				}
			if (length(jd) > 0){
				zigu2 <- c(zigu2, ulow[jd])
				zigv2 <- c(zigv2, vlow[jd])
				}
			}
		lonlatu <- t(A)%*%t(cbind(zigu1,zigv1))
		lonu <- lonlatu[1,]/60/c0 + lon0
		latu <- lonlatu[2,]/60 + lat0
		lonlatd <- t(A)%*%t(cbind(zigu2, zigv2))
		lond <- lonlatd[1,]/60/c0 + lon0
		latd <- lonlatd[2,]/60 + lat0
		ziglon1 <- lonu
		ziglat1 <- latu 
		ziglon2 <- lond
		ziglat2 <- latd
		xz1 <- (ziglon1-lon0)*60*cos(lat0/180*pi)
		yz1 <- (ziglat1-lat0)*60
		xz2 <- (ziglon2-lon0)*60*cos(lat0/180*pi)
		yz2 <- (ziglat2-lat0)*60
		dist1 <- sum(sqrt(diff(xz1)^2 + (diff(yz1))^2))
		dist2 <- sum(sqrt(diff(xz2)^2 + (diff(yz2))^2))
		out <- list(ziglon1 = ziglon1, ziglat1 = ziglat1, 
		ziglon2 = ziglon2, ziglat2 = ziglat2, 
		dist1 = dist1, dist2 = dist2, zigu1 = zigu1, zigv1 = zigv1)
		}