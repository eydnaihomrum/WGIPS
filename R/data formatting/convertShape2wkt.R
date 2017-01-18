install.packages("maptools")
install.packages("rgdal")
install.packages("rgeos")
library(rgeos)
library(maptools)
library(rgdal)
wd1 <- setwd("..\\Rekestrata\\") ## Link to data folder 
shp1 <- "shrimp_stratsys_NSSK_2015a" ## Shapefile (layer) 
out.wkt <- "shrimp_stratsys_NSSK_2015a.wkt"

xx <- readOGR(dsn=wd1,layer=shp1) 
plot(xx, border="blue", axes=TRUE, las=1) 
text(coordinates(xx), labels=row.names(xx), cex=0.6)

xx1 <- spTransform(xx,CRS("+proj=longlat +datum=WGS84"))
plot(xx1, border="red", axes=TRUE, las=1) 
names.strata <- as.character(xx1$Strat_desc)
x1 <- writeWKT(xx1, byid = TRUE)
x2 <- gsub("POLYGON","MULTIPOLYGON (",x1 )
x3 <- gsub("))",")))",x2)

for(i in 1:length(names.strata)){
  x3[i] <- sub("MULTIPOLYGON", paste(names.strata[i], "MULTIPOLYGON", sep="\t"), x3[i])
}

writeLines(x3,paste0(wd1,'\\',out.wkt))

