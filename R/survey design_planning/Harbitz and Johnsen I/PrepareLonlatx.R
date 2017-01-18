#### Prepare lonlatx.txt
# lonlatx: matrix with polygon node positions in decimal degrees
#  col.1: longitude, col.2: latitude, col.3: stratum number (no)
datapath <- "../Kartdata_2015/"
x1 <- read.table(paste0(datapath,"stratum1.txt"),header=T)
x1$Stratum <- 1
x1 <- x1[,c(2,1,3)]; x1[,1] <- x1[,1]

x2 <- read.table(paste0(datapath,"stratum2.txt"),header=T)
x2$Stratum <- 2
x2 <- x2[,c(2,1,3)]


x3 <- read.table(paste0(datapath,"stratum3.txt"),header=T)
x3$Stratum <- 3
x3 <- x3[,c(2,1,3)]

x4 <- read.table(paste0(datapath,"stratum4.txt"),header=T)
x4$Stratum <- 4
x4 <- x4[,c(2,1,3)]

x5 <- read.table(paste0(datapath,"stratum5.txt"),header=T)
x5$Stratum <- 5
x5 <- x5[,c(2,1,3)]

x.final <- rbind(x1,x2,x3,x4,x5)
write.table(x.final, file=paste0(datapath,"lonlatx.txt"), row.names = F, col.names=F)

# lonlat12x: matrix with one entrance point and one exit point in decimal degrees
# (the entrance and exit points defines head direction)
get.lonlat12x <- function(x){
  y1 <- y2 <- (min(x[,2])+max(x[,2]))/2
  x2 <- min(x[,1]) - 0.1
  x1 <- max(x[,1]) + 0.1
  out <- c(y1,x1,y2,x1)
}

p1 <- get.lonlat12x(x1)
p2 <- get.lonlat12x(x2)
p3 <- get.lonlat12x(x3)
p4 <- get.lonlat12x(x4)
p5 <- get.lonlat12x(x5)

ans1 <- as.data.frame(rbind(p1,p2,p3,p4,p5))
ans1$stratum <- 1:5
#write.table(ans1, file=paste0(datapath,"lonlat12x.txt"), row.names = F, col.names=F)


