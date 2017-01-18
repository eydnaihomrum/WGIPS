# zigzagrnd.r

# Calculates randomized zigzag transects for a given polygon and a given head direction
# The algorithm is based on equidistant lines across the head direction and the
# intersection between these lines and the polygon, with a random, uniformly distributed
# start line
#
# Reference:
# Strindberg, S. and Buckland, S.T., 2004. Zigzag Survey Designs in Line Transect Sampling.
# Journal of Agricultural, Biological, and Environmental Statistics. Vol. 9 No. 4 p443-461

# Script made by Alf Harbitz May 2010, e-post alf.harbitz@imr.no, wtel +4777609731

# NOTES:
# 1) The program may crash if the time t is too small, allowing only a very small number of transects
# 2) The calculation of polygon area with areapl requires libraries splancs and sp

# INPUT
# no = field number
# t = vessel time in hours
# vel = vessel velocity in knots
 
# OUTPUT
# out$cruiseline = a dataframe with:
#                 * lon = vector with decimal longitude values in degrees at polygon nodes
#                 * lat = vector with decimal latitudes values in degrees at polygon nodes  
# out$Area = Area of the selected Survey area
# out$Distance = Distance of the cruiseline
# out$Sur.cov = Degree of survey coverage (dist/sqrt(A))

# lon = vector with decimal longitude values in degrees at polygon nodes
# lat = vector with decimal latitudes values in degrees at polygon nodes
# out = mx2 matrix with random zigzag decimal longitude and latitude values in degrees
# lonrnd = 1st column of out
# latrnd = 2nd column of out

# FILES TO BE READ: 
# fieldname: field names with assigned field numbers
# lonlatx: matrix with polygon node positions in decimal degrees
#	col.1: longitude, col.2: latitude, col.3: field number (no)
# lonlat12x: matrix with one entrance point and one exit point in decimal degrees
# (the entrance and exit points defines head direction)
 
datapath <- "../Kartdata/"
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
tim <- readline('time available in hours,  t = ? ')
vel <- readline('vessel speed in knots,  vel = ? ')
nsim <- readline('no of runs,            nsim = ? ')
no <- as.numeric(no)
tim <- as.numeric(tim)
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
lb12 <- as.numeric(lb12)								# converts to numeric (lb12 is one line)
lon1 <- lb12[1]				# decimal longitude of entrance point
lat1 <- lb12[2]				# decimal latitude of entrance point
lon2 <- lb12[3]				# decimal longitude of exit point
lat2 <- lb12[4]				# decimal latitude of exit point
n1 <- length(lon);			# number of nodes in polygon 
lon <- matrix(lon,n1,1)		# ascertains lon to be column vector
lat <- matrix(lat,n1,1)		# ascertains lat to be column vector
 
# ascertains lon(1) = lon(end) and lat(1) = lat(end) if the difference of entrance and exit points
# exceed 1e-6 deg for either lon-values, lat-values or both
if (abs(lon[1]-lon[n1]) > 1e-6 | abs(lat[1]-lat[n1]) > 1e-6)	# TRUE if difference exceeds threshold
	{
    lon <- c(lon, lon[1]);	# adds a new last lon-value identical to first one
    lat <- c(lat, lat[1]);	# adds a new last lat-value identical to first one
    n1 <- n1 + 1;				# adjusts n to define the length of lon and lat
    }
lon <- matrix(lon,n1,1)		# ascertains lon to be a vertical vector
lat <- matrix(lat,n1,1)		# ascertains lat to be a vertical vector



## Espen has removed this if
#if (nsim == 1)
#	{
#	plot(lon,lat,type = 'l')# plots polygon
#	}
## ===============
## Espens line
plot(lon,lat,type = 'l',main=name.field, col=2)
#points(c(lon1,lon2),c(lat1,lat2),col=2)
text(c(lon1,lon2),c(lat1,lat2),labels=1:2,col=2)


# TRANSFORMATION FROM (lon,lat) TO CARTESIAN NAUTICAL MILE COORDINATES (x,y)
fac <- pi/180					# trigonometric factor from deg to rad
lonm <- mean(lon)			# mean of polygon longitude points
latm <- mean(lat)			# mean of polygon latitude points
c0 <- cos(latm*fac) 		# factor to get cartesian x-coordinates from longitude
x <- (lon-lonm)*c0*60		# polygon x-values in nautical miles 
y <- (lat-latm)*60			# polygon y-values in nautical miles
x1 <- (lon1-lonm)*c0*60		# (x1,y1) = entrance point in cartesian coordinates
y1 <- (lat1-latm)*60         
x2 <- (lon2-lonm)*c0*60		# (x2,y2) = exit point in cartesian coordinates
y2 <- (lat2-latm)*60

# TRANSFORMATION (ROTATION) FROM (x,y) TO (u,v) WITH u = HEAD DIRECTION AXIS
tet <- atan2(y2-y1,x2-x1)	# head direction
# rot = rotation matrix
rot <- matrix(cbind(cos(tet), sin(tet),-sin(tet), cos(tet)),2,2, byrow = TRUE)
uv <- rot%*%matrix(t(cbind(x,y)),2,n1)	# 2xn matrix with u and v values in 1st and 2nd row
u1v1 <- rot%*%matrix(cbind(x1, y1, x2, y2),2,2,byrow=TRUE)
u <- uv[1,1:n1];				# (u,v) = polygon nodes
v <- uv[2,1:n1];
u <- matrix(u,n1,1)
v <- matrix(v,n1,1)
u1 <- u1v1[1,1]				# (u1,v1) = entrance point
v1 <- u1v1[1,2]
u2 <- u1v1[2,1]				# (u2,v2) = exit point
v2 <- u1v1[2,2]
browser()
# ascertains no vertical polygon segments in (u,v)
j0 <- diff(u)*diff(v)
if (length(which(j0 == 0)) > 0)			# TRUE if at least one vertical segment
	{
    u <- u + 1e-6*sd(u)*(runif(n1)-.5);	# slight adjustment of all polygon nodes
    v <- v + 1e-6*sd(v)*(runif(n1)-.5);    
    u <- vector(u,n1,1);
    v <- vector(v,n1,1)
    }
     
# CALCULATES CONSTANT u-DISTANCE c BETWEEN SUCCEEDING ZIGZAG TURNING POINTS
poly <- cbind(x,y)	# makes a nx2 matrix with x and y in 1st and 2nd columns
A <- areapl(poly)	# calculates area of polygon in square nautical miles
umin <- min(u)		# minimum u-value of polygon
umax <- max(u)		# maximum u-value of polygon
c <- A/sqrt((vel*tim)^2 - (umax-umin)^2)  # u-spacing between zig-zag turning points

# CALCULATES LINE COEFFICIENTS (1,b), v = a + bu, FOR EACH POLYGON SEGMENT
uf <- u[1:n1-1]		# u-values of startpoint of polygon segments
ul <- u[2:n1]			# u-values of endpoint of polygon segments
vf <- v[1:n1-1]		# v-values of startpoint of polygon segments
vl <- v[2:n1]			# v-values of endpoint of polygon segments
du <- ul - uf			# difference in u between end and start of segments
dv <- vl - vf			# difference in v between end and start of segments
b <- dv/du	      	# slope of polygon segments
a <- vl - b*ul	  	# intercept of polygon segments

for (jsim in 1:max(1,nsim)) # replicate loop
{
# CALCULATE ZIGZAG u-VALUES urnd WITH RANDOM START POINT, 
# min(urnd) < umin < umax < max(urnd)
# ONLY ONE urnd-POINT AT EACH SIDE OF THE POLYGON EXTENSION IN u-DIRECTION 
urnd <- umin + 2*runif(1)*c	# random reference u-point for zig zag points
urnd <- seq(urnd,umax,by=c)	# adding u-values for turning points until umax
urnd <- c(urnd, max(urnd)+c)	# adding one u-value to the right of polygon
urnd <- c(urnd[1]-c, urnd)		# adding one u-value to the left
if (urnd[1] > umin)
	{	
    urnd <- c(urnd[1]-c, urnd)	# ascertains one u-point left of polygon
    }
nurnd <- length(urnd)           	# length of u-vector with turning points

# CALCULATES UPPER AND LOWER v-VALUES FOR INTERSECTION BETWEEN ZIGZAG TRANSECTS
# AND POLYGON
vlow <- matrix(0,nurnd,1)		# initiating vector with lower v-values at turning points
vupp <- vlow						# initiating vector with upper v-values at turning points
for (j in 2:(nurnd-1))			# only nurnd - 2 values inside polygon
	{
    uj <- urnd[j]				# j'th u-value
    duj <- (uj-uf)*(uj-ul)		# criterion for finding appropriate polygon segments
    k <- which(duj < 0)			# index to segments including uj
    vk <- a[k] + b[k]*uj		# v-values at uj for segments including uj
    vlow[j] <- min(vk)			# assigning the lowest v-value to vlow
    vupp[j] <- max(vk)			# assigning the largest v-value to vupp
    }
vlow[1] <- vlow[2]				# assigning the "neighbour" value at first u-value
vupp[1] <- vupp[2]              
vlow[nurnd] <- vlow[nurnd-1]	# assigning the "neighbour" value at last u-value
vupp[nurnd] <- vupp[nurnd-1]

# CALCULATES ONE ZIGZAG WITH LOWER START AND ONE WITH UPPER START (LEFT) 
j1 <- seq(1,nurnd,2)			# index at either upper or lower turning points
ind <- 1:nurnd					# index vector
j2 <- setdiff(ind,j1)			# index at turning points not identified by j1
vrnd1 <- matrix(0,nurnd,1)		# v-values including turning points when start low-left
vrnd2 <- matrix(0,nurnd,1)		# v-values including turning points when start upper-left
vrnd1[j1] <- vlow[j1]			# low v-values when starting low
vrnd1[j2] <- vupp[j2]			# large v-values when starting low
vrnd2[j1] <- vupp[j1]			# large v-values when starting high
vrnd2[j2] <- vlow[j2]			# low v-values when staring high
urnd0 <- urnd						# maintains reference urnd-vector

# FINDS ENTRANCE AND EXIT INTERSECTION BETWEEN ZIGZAG AND POLYGON  

# Intersection between lower start zigzag and first polygon segment
vrnd <- vrnd1		# vrnd = help variable
errfirst1 <- 0	# errfirst = 1 if no intersection with first polygon segment
b1 <- (vrnd[2]-vrnd[1])/(urnd0[2]-urnd0[1])		# slope of first zig-zag transect
a1 <- vrnd[2] - b1*urnd0[2]						# intersect of first zig-zag transect
uu <- -(a1 - a)/(b1 - b)							# u-value at intersection with polygon segments
ind <- which((uu - ul)*(uu - uf) < 0)			# index to which segments can be intersected
vv <- a[ind] + b[ind]*uu[ind]						# actual vv-values (2)
vv1 <- vv[1]											# Assumes first vv-value is correct
uu1 <- uu[ind[1]]									# Assumes first uu-value is correct
if (abs(vv1 - vrnd[2]) < abs(vv[2] - vrnd[2]))	# criterion for second value
    {
    vv1 <- vv[2];				# Assigns second value
    uu1 <- uu[ind[2]]
    }
if (uu1 > urnd0[2])
	{
    errfirst1 <- 1
    }
urnd[1] <- uu1
vrnd[1] <- vv1
vrnd1 <- vrnd
 
# Intersection between lower start zigzag and last polygon segment
vrnd <- vrnd1		# vrnd = help variable adjusted for possible intersection with first segment
errlast1 <- 0		# errlast = 1 if no intersection with last polygon segment
b1 <- (vrnd[nurnd]-vrnd[nurnd-1])/(urnd0[nurnd]-urnd0[nurnd-1]) # slope of last zig-zag transect
a1 <- vrnd[nurnd] - b1*urnd0[nurnd]		# intersect of last zig-zag transect
uu <- -(a1 - a)/(b1 - b)					# u-value at intersection with polygon segments
ind <- which((uu - ul)*(uu - uf) < 0)	# index to which segments can be intersected    
vv <- a[ind] + b[ind]*uu[ind]				# actual vv-values (2)
vv1 <- vv[1]									# Assumes first vv-value is corrcet
uu1 <- uu[ind[1]]							# Assumes first uu-value is correct
if (abs(vv1 - vrnd[nurnd-1]) < abs(vv[2] - vrnd[nurnd-1])) # criterion for second value
    {
    vv1 <- vv[2]; 
    uu1 <- uu[ind[2]]		   # Assigns second value
    }
if (uu1 < urnd0[nurnd-1])
    {
    errlast1 <- 1
    }
urnd[nurnd] <- uu1;
vrnd[nurnd] <- vv1;
vrnd1 <- vrnd
urnd1 <- urnd

# Intersection between upper start zigzag and first polygon segment
vrnd <- vrnd2		# vrnd = help variable
errfirst2 <- 0	# errfirst2 = 1 if no intersection with first polygon segment
b1 <- (vrnd[2]-vrnd[1])/(urnd0[2]-urnd0[1]) # slope of first zig-zag transect
a1 <- vrnd[2] - b1*urnd0[2]				# intersect of first zig-zag transect
uu <- -(a1 - a)/(b1 - b)					# u-value at intersection with polygon segments
ind <- which((uu - ul)*(uu - uf) < 0)	# index to which segments can be intersected
vv <- a[ind] + b[ind]*uu[ind]				# actual vv-values (2)
vv1 <- vv[1]									# Assumes first vv-value is correct
uu1 <- uu[ind[1]]							# Assumes first uu-value is correct
if (abs(vv1 - vrnd[2]) < abs(vv[2] - vrnd[2]))    # criterion for second value
    {
    vv1 <- vv[2];	        # Assigns second value
    uu1 <- uu[ind[2]]
    }
if (uu1 > urnd0[2]) 
	{
    errfirst2 <- 1
    }
urnd[1] <- uu1
vrnd[1] <- vv1
vrnd2 <- vrnd
vrnd <- vrnd2
 
# Intersection between upper start zigzag and last polygon segment
vrnd <- vrnd2		# vrnd = help variable adjusted for possible intersection with first segment
errlast2 <- 0		# errlast2 = 1 if no intersection with last polygon segment
b1 <- (vrnd[nurnd]-vrnd[nurnd-1])/(urnd0[nurnd]-urnd0[nurnd-1]) # slope of first zig-zag transect
a1 <- vrnd[nurnd] - b1*urnd0[nurnd]		# intersect of last zig-zag transect
uu <- -(a1 - a)/(b1 - b)					# u-value at intersection with polygon segments
ind <- which((uu - ul)*(uu - uf) < 0)	# index to which segments can be intersected   
vv <- a[ind] + b[ind]*uu[ind]				# actual vv-values (2)
vv1 <- vv[1]									# Assumes first vv-value is corrcet
uu1 <- uu[ind[1]]							# Assumes first uu-value is correct
if (abs(vv1 - vrnd[nurnd-1]) < abs(vv[2] - vrnd[nurnd-1])) # criterion for second value
    {
    vv1 <- vv[2];	        # Assigns second value
    uu1 <- uu[ind[2]]
    }
if (uu1 < urnd0[nurnd-1])
	{
    errlast2 <- 1
    }
urnd[nurnd] <- uu1
vrnd[nurnd] <- vv1
vrnd2 <- vrnd;
urnd2 <- urnd;

# adjusts the zig zag for entrance and exit points
if (errfirst1 > 0) 
    {
    urnd1 <- urnd1[2:nurnd]; 
    vrnd1 <- vrnd1[2:nurnd]
    }
if (errlast1 > 0)
	{
    urnd1 <- urnd1[1:length(urnd1)-1]; 
    vrnd1 <- vrnd1[1:length(vrnd1)-1]
    }
nurnd1 <- length(urnd1); 
if (errfirst2 > 0)
	{
    urnd2 <- urnd2[2:nurnd]; 
    vrnd2 <- vrnd2[2:nurnd]
    }
if (errlast2 > 0) 
	{
    urnd2 <- urnd2[1:length(urnd2)-1]; 
    vrnd2 <- vrnd2[1:length(vrnd2)-1]
    }
nurnd1 <- length(urnd1)
nurnd2 <- length(urnd2);
urnd1 <- matrix(urnd1,nurnd1,1)
vrnd1 <- matrix(vrnd1,nurnd1,1)
urnd2 <- matrix(urnd2,nurnd2,1)
vrnd2 <- matrix(vrnd2,nurnd2,1)
 
# BACKTRANSFORMATION FROM (u,v) TO (x,y) AND (longitude,latitude)
xy1 <- t(rot)%*%matrix(t(cbind(urnd1,vrnd1)), 2, nurnd1)		# start low
xy2 <- t(rot)%*%matrix(t(cbind(urnd2,vrnd2)),2, nurnd2)		# start high
x1 <- xy1[1,]	      # (x1,y1) = zig-zag nodes with start low
y1 <- xy1[2,]
x2 <- xy2[1,]	      # (x2,y2) = zig-zag nodes with start high
y2 <- xy2[2,]
hrlow <- sum(sqrt(diff(x1)^2 + diff(y1)^2))/vel	# actual t by zigzag with low start
hrupp <- sum(sqrt(diff(x2)^2 + diff(y2)^2))/vel	# actual t by zigzag with high start
lonrnd <- x1/c0/60 + lonm     	# (lonrnd,latrnd) = zig-zag nodes (longitude,
latrnd <- y1/60 + latm        	# latitude) with start low
lonrnd1 <- lonrnd
latrnd1 <- latrnd
lonrnd2 <- x2/c0/60 + lonm      	# (lonrnd,latrnd) = zig-zag turning points (longitude,
latrnd2 <- y2/60 + latm         	# latitude) with start high
ru <- runif(1)	               	# random uniform U(0,1) value

# PICKS RANDOM ZIGZAG TRANSECT BETWEEN THW LOW START AND HIGH START ALTERNATIVES
if (ru > .5)                  	# picks randomly if start low or high is chosen
	{
    lonrnd <- lonrnd2;         	# ru < .5: start low, ru > .5: start high (from left)
    latrnd <- latrnd2
    }
lines(lonrnd,latrnd)
}

# output: nx2 matrix with random zigzag decimal longitudes and latitudes
dist1 <- 0
for(i in 1:(length(lonrnd)-1)){
  dist1 <- dist1 + (spDistsN1(cbind(lonrnd,latrnd),cbind(lonrnd,latrnd)[i,], longlat=T)[i+1])
  }
dist1 <- dist1/1.852  

#dist1 <- spDistsN1(cbind(lonrnd,latrnd),cbind(lonrnd,latrnd)[1,], longlat=T)[i+1]
browser()
polyg <- cbind(lon,lat)
out <- list(name.field=name.field,
        cruiseline = as.data.frame(cbind(lonrnd,latrnd)), 
        area.nm2 <- areaPolygon(polyg)/1852/1852, # Meters2,
        Distance = dist1,
        Sur.cov = dist1/ sqrt(A)) 
        
        




## Create file to MaxSea
## Create ascii fil til import i MaxSea (se manual for code)
item.type <- 257  # Linje
item.id <- 20      # Setter id-coder på linjene
item.col <- 1     # Blå
north <- "N"
west <- "E"

##
filnavn.ut <- paste("../Kartdata/rutelinjer/",name.field,".asc",sep="")
ant.dec.pos <- 3
n.row <- nrow(out$cruiseline)
                      
out1 <- cbind(rep(item.type,n.row),
            rep(item.id,n.row),
            rep(item.col,n.row),
            round(out$cruiseline$latrnd,ant.dec.pos),
            rep(north,n.row),
            round(out$cruiseline$lonrnd,ant.dec.pos),
            rep(west,n.row))
            
write.table(out1,file=filnavn.ut,row.names=F, col.names=F,sep=",",quote=F)

polyg <- cbind(lon,lat)
polyg.sp = SpatialPolygons(list(Polygons(list(Polygon(polyg)), "x")))
#plot(polyg.sp)
#points(spsample(polyg.sp, n = 10, "random"), pch = 2, col=2)
spsample(polyg.sp, n = 10, "random")

  tmp <- as.data.frame(cbind(floor(out$cr$lonrnd), round(60*(out$cr$lonrnd- floor(out$cr$lonrnd)),2), 
    floor(out$cr$latrnd), round(60*(out$cr$latrnd- floor(out$cr$latrnd)),2)))
  names(tmp) <- c("LonDeg","LonMin","LatDeg","LatMin")  
  name.field
  vel # Speed
  dist1 <- out$Distance## Sailing distance
  area.nm2 <- areaPolygon(polyg)/1852/1852 # Meters2
  Sur.cov <- dist1/ sqrt(area.nm2)
  print(tmp)  
  ## TABLE
  filnavn1 <- paste("../Kartdata/rutelinjer/INFO",name.field,".txt",sep="")
  capture.output( cat("\n",format(c(date()),width=20, justify = "left")), file=filnavn1)
  capture.output( cat("\n", "Stratum (Toktområde)    ", format(name.field,width=7, justify ="right")), file=filnavn1, append=T)
  capture.output( cat("\n", "Speed and time available", format(c(paste(vel,"knots"),paste(tim, "h")),width=7, justify ="right")), file=filnavn1, append=T)
  capture.output( cat("\n", "Stratum area (n.mi2)    ", format(area.nm2,width=7, justify ="right")), file=filnavn1, append=T)
  capture.output( cat("\n", "Sailing distance (n.mi) ", format(dist1,width=7, justify ="right")), file=filnavn1, append=T)
  capture.output( cat("\n", "Survey coverage         ", format(Sur.cov,width=7, justify ="right")), file=filnavn1, append=T)
  capture.output( cat("\n", " "), file=filnavn1, append=T)
  capture.output( cat("\n", "Transect positions      "), file=filnavn1, append=T)
  capture.output( cat("\n", " "), file=filnavn1, append=T)
  capture.output(cat("\n", " "), tmp, file= filnavn1, append=T)  

  filnavn2 <- paste("../Kartdata/rutelinjer/TRACK",name.field,".txt",sep="")
  out2 <- as.data.frame(cbind(1,out$cr))
  names(out2) <- c("Line","Longitude", "Latitude")
  write.csv(out2, file=filnavn2,row.names=F)