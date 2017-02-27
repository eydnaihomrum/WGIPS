
##################################
## Datamining - quality checklist before uploading data to pgnapes or ices-db 
## 
## Example data: Faroese logbook-, catch- and biology IBWSS data 2011 were extracted from pgnapes-access-database
## OBS! rename columns and do other alterations according to national format
##
## author: Eydna í Homrum, Faroe Marine Research Institute, eydnap@hav.fo
##################################


rm(list = ls())

#### Packages ####
library(mapplots) # draw.pie
data("coast") # coastlines
library(maps)
library(mapdata)
library(gdata)  # drop.levels
library(reshape) # melt and cast


#### import data #### 
## read tables
logbook <- read.table("d:/wgips/2017/checklist/mh1111/mh1111_logbook.txt", sep = ",", header = TRUE)
names(logbook) <- tolower(names(logbook))
str(logbook)

biology <- read.table("d:/wgips/2017/checklist/mh1111/mh1111_biology.txt", sep = ",", header = TRUE)
names(biology) <- tolower(names(biology))
str(biology)

### OBS!! Erroneous species in example dataset - DLF should be TPA
biology$species <- as.character(biology$species)
biology$species[biology$species == 'DLF'] <- 'TPA'
biology$species <- as.factor(biology$species)

# reference data: 
# here Faroese IBWSS biology data for last 2 years (= 2009 and 2010)
# alternatively - query data based on a lon,lat polygon in the relevant area and season for 2,5 or 10 years ....

two <- read.table("d:/wgips/2017/checklist/mh1111/ibwss_biology_fo2ar.csv", sep = ",", header = TRUE)
names(two) <- tolower(names(two))
str(two)

catch <- read.table("d:/wgips/2017/checklist/mh1111/mh1111_catch.txt", sep = ",", header = TRUE)
names(catch) <- tolower(names(catch))
str(catch)

acoustic <- read.table("d:/wgips/2017/checklist/mh1111/mh1111_acoustic.txt", sep = ",", header = TRUE)
names(acoustic) <- tolower(names(acoustic))
str(acoustic)

acousticvalues <- read.table("d:/wgips/2017/checklist/mh1111/mh1111_acousticvalues.txt", sep = ",", header = TRUE)
names(acousticvalues) <- tolower(names(acousticvalues))
str(acousticvalues)

#if applicable
plankton <- read.table("d:/wgips/2017/checklist/mh1111/mh1111_plankton.txt", sep = ",", header = TRUE)
names(plankton) <- tolower(names(plankton))
str(plankton)

# species table - to convert FAO codes to latin names
spec <- read.table("d:/wgips/2017/checklist/pgnapes_talvur/pgnapes_species.csv", sep = ",", header = TRUE)
names(spec) <- tolower(names(spec))
str(spec)
## These are the species in the pgnapes-database
## if you encounter other species - please contact the database responsible person


#### Step 1: Biology table ####
## format - one record per fish

# fishlength 'length' - fishweight 'weight'

# if fishlength is in .1 cm or half cm - create column with 1cm-length group 'lgroup'
biology$lgroup <- floor(biology$length)
# Fulton K - from fishlength 'length' and fishweight 'weight'
biology$k <- biology$weight/biology$length^3*100
two$k <- two$weight/two$length^3*100


# Age is = ageotolith - if agescale is filled in (e.g. Herring NO and IS) agescale is used 
biology$age <- biology$ageotolith
biology$age[!is.na(biology$agescale)] <- biology$agescale[!is.na(biology$agescale)]
table(biology$species[!is.na(biology$age)])

# Age for reference data 
two$age <- two$ageotolith
two$age[!is.na(two$agescale)] <- two$agescale[!is.na(two$agescale)]
table(two$species[!is.na(two$age)])

#### Step 1.1 Length frequencies ####

# set mfrow according to number of species
length(unique(biology$species))

#png("d:/wgips/2017/checklist/lengthdistribution_test.png", width = 16, height = 12, units = "cm", res = 300, pointsize = 10)

# mfrow = c(3,3) - change if more than 9 species with lengths ...
opar <- par(mfrow = c(3,3), mar = c(2,2,1,1), oma = c(2,2,3,0))
for(i in sort(unique(biology$species))){
    dd <- biology[biology$species == i,]
    (ldist <- aggregate(x = dd$lgroup, by = list(dd$lgroup), FUN = length))
    (ldist2 <- merge(x = data.frame(length = seq(min(dd$lgroup)-3, max(dd$lgroup)+3,1)), y = ldist, by.x = "length", by.y = "Group.1", all.x = TRUE))
    (ldist2$procent <- ldist2$x/sum(ldist2$x, na.rm = TRUE)*100)
    
    mp <- barplot(ldist2$procent, ylim = c(0, 1.1*max(ldist2$procent, na.rm = TRUE)), las = 2, xpd = FALSE)
    axis(side = 1, labels = seq(min(dd$lgroup)-3, max(dd$lgroup)+3,1), at = mp)
    box()
    legend("topleft", legend = paste(as.character(i), '\n', spec$scientific_name[spec$speciesid == as.character(i)], sep = ""), bty = "n") 
    legend("topright", legend = paste("n = ", as.character(sum(ldist2$x, na.rm = TRUE))), bty = "n") 
}
mtext(side = 3, outer = TRUE, "Length distribution")
mtext(side = 2, outer = TRUE, "Frequency (%)", line = .75)

#dev.off()

#### Step 1.2 - length-weight relationship ####



#png("d:/wgips/2017/checklist/weight_vs_length_test.png", width = 16, height = 12, units = "cm", res = 300, pointsize = 10)

# mfrow = c(3,4) - change if more than 12 species with lengths ...
opar <- par(mfrow = c(3,3), mar = c(2,2,.1,2), oma = c(0,0,3,0))
for(i in sort(unique(biology$species))){
  aa <- biology[biology$species == i,]
  bb <- two[two$species == i,]
  (xmax <- c(max(aa$length), max(bb$length)))
  (ymax <- c(max(aa$weight, na.rm = TRUE), max(bb$weight, na.rm = TRUE)))
  plot(bb$weight ~ bb$length, axes = FALSE, pch = 3, col = "grey",
       xlim = c(0, max(xmax)), ylim = c(0, max(ymax)))
  
  points(aa$weight ~ aa$length, pch = 4, col = "red")
  axis(side = 1)
  axis(side = 2)
  box()
  legend("topleft", paste(as.character(i), '\n', spec$scientific_name[spec$speciesid == as.character(i)], sep = ""), bty = "n")
}
mtext(side = 3, outer = TRUE, "Length-weight relationship", line = 1.5)
par(opar)
#dev.off()

## If some species only have FAO three letter code - and are missing the latin name, then this speciesID is not in the species-table

## overview over species with weight information 
## - obs! no weight measurements for 'LPR' and 'PLS'
table(biology$species[!is.na(biology$weight)])

## Please check:
#  length has to be in 'cm'!
#  weight has to be in 'g'!

#### Step 1.3 - length-condition factor relationship ####

#png("d:/wgips/2017/checklist/condition_vs_length_test.png", width = 16, height = 12, units = "cm", res = 300, pointsize = 10)

opar <- par(mfrow = c(3,3), mar = c(2,2,.1,2), oma = c(0,0,3,0))
for(i in sort(unique(biology$species))){
  cc <- biology[biology$species == i,]
  dd <- two[two$species == i,]
  (xmax <- c(max(cc$length), max(dd$length)))
  (ymax <- c(max(cc$k, na.rm = TRUE), max(dd$k, na.rm = TRUE)))
  plot(dd$k ~ dd$length, axes = FALSE, pch = 3, col = "grey", 
       xlim = c(0, max(xmax)), ylim = c(0, max(ymax)))
  points(cc$k ~ cc$length, pch = 4, col = "red")
  axis(side = 1)
  axis(side = 2)
  box()
  legend("topleft", paste(as.character(i), '\n', spec$scientific_name[spec$speciesid == as.character(i)], sep = ""), bty = "n")
}
mtext(side = 3, outer = TRUE, "Length-condition factor relationship", line = 1.5)
par(opar)
#dev.off()

#### Step 1.4 - length - age relationship ####

#png("d:/wgips/2017/checklist/condition_vs_length_test.png", width = 16, height = 12, units = "cm", res = 300, pointsize = 10)

# mfrow = c(2,2) - change if more than 4 species with ages
opar <- par(mfrow = c(2,2), mar = c(2,2,.1,2), oma = c(0,0,3,0))
for(i in sort(unique(biology$species[!is.na(biology$age)]))){
  ee <- biology[biology$species == i & !is.na(biology$age),]
  ff <- two[two$species == i & !is.na(two$age),]
  plot(ff$length ~ ff$age, axes = FALSE, pch = 3, col = "grey", 
       xlim = c(0, max(max(ee$age), max(ff$age))), ylim = c(0, max(max(ee$length), max(ff$length))))
  points(ee$age, ee$length, col = "red", pch = 4)
  axis(side = 1)
  axis(side = 2)
  box()
  legend("bottomright", paste(as.character(i), '\n', spec$scientific_name[spec$speciesid == as.character(i)], sep = ""), bty = "n")
}
mtext(side = 3, outer = TRUE, "Length-age relationship", line = 1.5)
par(opar)
#dev.off()

#### Step 1.5 - Age distribution  ####

#png("d:/wgips/2017/checklist/agedist_test.png", width = 16, height = 12, units = "cm", res = 300, pointsize = 10)

# mfrow = c(2,2) - change if more than 4 species with ages ...
opar <- par(mfrow = c(2,2), mar = c(2,2,1,1), oma = c(2,2,3,0))
for(i in sort(unique(biology$species[!is.na(biology$age)]))){
    dd <- biology[biology$species == i,]
    (adist <- aggregate(x = dd$age, by = list(dd$age), FUN = length))
    (adist2 <- merge(x = data.frame(age = 1:max(dd$age, na.rm = TRUE)), y = adist, by.x = "age", by.y = "Group.1", all.x = TRUE))
    (adist2$procent <- adist2$x/sum(adist2$x, na.rm = TRUE)*100)
    
    mp <- barplot(adist2$procent, ylim = c(0, 1.1*max(adist2$procent, na.rm = TRUE)), las = 2)
    axis(side = 1, labels = 1:max(dd$age, na.rm = TRUE), at = mp)
    box()
    legend("topleft", legend = paste(as.character(i), '\n', spec$scientific_name[spec$speciesid == as.character(i)], sep = ""), bty = "n") 
    legend("topright", legend = paste("n = ", as.character(sum(adist2$x, na.rm = TRUE))), bty = "n") 
}
mtext(side = 3, outer = TRUE, "Age distribution")
mtext(side = 2, outer = TRUE, "Frequency (%)", line = .75)
par(opar)

#dev.off()


#### Step 2 - Catch table ####

#### Step 2.1 - Species composition and total catch ####
str(logbook)
str(catch)

## main species on map 
# 1 - combine catch data with logbook data
catpos <- merge(catch, logbook[,c('station', 'lon', 'lat')], by.x = 'station', by.y = 'station', all = TRUE)
str(catpos)
# 2 - subset to species to plot
mspec <- drop.levels(catpos[catpos$species %in% c('HER', 'MAC', 'WHB'),])
str(mspec)

# 3 - create lists to plot - from package 'mapplots'
hh <- make.xyz(x = mspec$lon, y = mspec$lat, z = mspec$catch, group = mspec$species)
str(hh)

opar <- par(mfrow = c(1,1), mar = c(2,2,1,1), oma = c(2,2,3,0))
map('worldHires', xlim = c(min(hh$x)-2, max(hh$x)+1), ylim = c(min(hh$y)-.5, max(hh$y)+.5))
box()
axis(side = 2, las = 2)
draw.pie(x = hh$x, y = hh$y, z = hh$z, radius = 0.5, col = rainbow(2))
legend.pie(min(hh$x)-1, max(hh$y)+.25, labels = dimnames(hh$z)[[2]], col = rainbow(2), radius = 0.2, bty = "n")
par(opar)

## Please check:
# Catch has to be in kg!

range(catch$catch)
#overview - catch by station and species
tapply(catch$catch, list(catch$station, catch$species), sum)

#### Step 2.2 - check station information in catch-table ####
## Towtime, wirelength, towspeed and trawldepth
##  (in pgnapes - is this also relevant for ices-db?? )
aggregate(x = catch$station, by = list(station = catch$station, sttype = catch$sttype, towtime = catch$towtime, wirelength = catch$wirelength, towspeed = catch$towspeed, trawldepth = catch$trawldepth, station = catch$station), FUN = length)

#### step 3 - Logbook table ####
str(logbook)

## station type 
table(logbook$station, logbook$sttype)
aggregate(logbook$station, by = list(sttype = logbook$sttype, station = logbook$station), FUN = length)

## wind direction
aggregate(logbook$station, by = list(windir = logbook$windir, station = logbook$station), FUN = length)

## wind speed
aggregate(logbook$station, by = list(winspeed = logbook$winspeed, station = logbook$station), FUN = length)


#### step 4 (and 5) - Acoustic and AcousticValue tables ####

str(acoustic)
str(acousticvalues)

acvalpos <- merge(x = acoustic, y = acousticvalues, by = 'log')
acval <- aggregate(x = list(sa = acvalpos$sa), by = list(log = acvalpos$log, aclon = acvalpos$aclon, aclat = acvalpos$aclat, species = acvalpos$species), FUN = sum)

range(acval$x[acval$species == 'HER'])
range(acval$x[acval$species == 'WHB'])

## logbook stations on top of cruisetrack
opar <- par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(0,0,2,0))
map('worldHires', col = "white", xlim = c(min(acoustic$aclon)-1, max(acoustic$aclon)+1), ylim = c(min(acoustic$aclat)-.5, max(acoustic$aclat)+.5))
box()
axis(side = 2, las = 2)
points(logbook$lon, logbook$lat, cex = 1.2, col = "darkgrey", pch = 16, type = "b")
points(acoustic$aclon, acoustic$aclat, pch = 46, col = "red")
draw.shape(coast, col = "darkgreen")
legend("topleft", legend = "Logbook stations on cruisetrack", bty = "n")
par(opar)

## sA-values on cruisetrack
opar <- par(mfrow = c(1,1), mar = c(1,1,1,1), oma = c(0,0,2,0))
for(i in c('WHB')){ ## set relevant species here e.g. c('WHB', 'HER',...)
  jj <- acval[acval$species == i,]
  
  map('worldHires', col = "white", xlim = c(min(acoustic$aclon)-1, max(acoustic$aclon)+1), ylim = c(min(acoustic$aclat)-.5, max(acoustic$aclat)+.5))
  box()
  axis(side = 2, las = 2)
  points(jj$aclon, jj$aclat, cex = log(jj$sa)/2, col = "darkgrey")
  points(acoustic$aclon, acoustic$aclat, pch = 46, col = "red")
  draw.shape(coast, col = "darkgreen")
  legend("topleft", legend = as.character(i), bty = "n")
}

mtext(side = 3, outer = TRUE, "Sum of sa by log-interval")
par(opar)


#### step 6 - Plankton table ####

str(plankton)
# 1 - combine planktondata and logbook data
plkpos <- merge(plankton, logbook[,c('station', 'lon', 'lat')], by.x = 'station', by.y = 'station', all = TRUE)
str(plkpos)

#### Step 6.1 Sumdryweight ####

opar <- par(mfrow = c(1,1))
map('worldHires', xlim = c(min(plkpos$lon)-1, max(plkpos$lon)+1), ylim = c(min(plkpos$lat)-.5, max(plkpos$lat)+.5))
box()
axis(side = 2, las = 2)
points(plkpos$lon, plkpos$lat, cex = log(plkpos$sumdrywt)/2, pch = 16, col = "red")
legend("topright", legend =  c(min(plkpos$sumdrywt, na.rm = TRUE), max(plkpos$sumdrywt, na.rm = TRUE)), pt.cex = log(c(min(plkpos$sumdrywt, na.rm = TRUE), max(plkpos$sumdrywt, na.rm = TRUE)))/2, pch = 16, col = "red")
par(opar)

#### Step 6.2 Sizefractioned zooplankton ####

## if zooplankton has been SIZEFRACTIONED, create file with zooplankton groups as variable 
## to draw pieplot of plankton groups
## OBS!!! - not applicable for Faroese data - not certain this is meaningful - feedback from countries that sizefraction zooplankton??

# melt from package 'reshape' - reconstruct data from several columns into one
plkmelt <- melt(plankton, id = c('station', 'country', 'vessel', 'cruise', 'sttype', 'year', 'uppstatdepth', 'lowstatdepth'), measured = c('sumdrywt', 'frac2000', 'frac1000', 'frac180', 'krill', 'fish', 'shrimp'))

plkmeltpos <- merge(plkmelt, logbook[logbook$sttype == 'WP2',c('station', 'lon', 'lat')], by.x = 'station', by.y = 'station', all = TRUE)
str(plkmeltpos)

(ii <- drop.levels(plkmeltpos[!(plkmeltpos$variable == 'sumdrywt'),]))
(ii <- ii[!is.na(ii$station),])
(iii <- make.xyz(ii$lon, ii$lat, ii$value, ii$variable))
str(iii)

opar <- par(mfrow = c(1,1))
map('worldHires', xlim = c(min(iii$x)-1, max(iii$x)+1), ylim = c(min(iii$y)-.5, max(iii$y)+.5))
box()
axis(side = 2, las = 2)
draw.pie(x = iii$x, y = iii$y, z = iii$z, radius = 0.25, col = rainbow(6))
legend.pie(min(iii$x), min(iii$y)-.1, labels = dimnames(iii$z)[[2]], col = rainbow(6), radius = 0.2, bty = "n")
par(opar)

## Please check:
## - dryweight has to be in mg!

range(plankton$sumdrywt, na.rm = TRUE)
## Probably very institute-specific, how data are stored - 
## but: - values must be related to area of net-opening (.25 m2) (mg dryweight/m2)

