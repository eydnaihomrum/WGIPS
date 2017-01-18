####
x1000 <- read.table("../Manifold/Depth_1000m.txt", header=T)
x20 <- read.table("../Manifold/Depth_20m.txt", header=T)

x11()
plot(x20)
points(x1000, col=1)
#### Stratum 1, grense i nord til EEZ, i sør til 7 15' (7.25)
x1 <- x1000[1677:1841,] ##Grense
x1b <- x20[1938:2248,]
s1 <- rbind(x1,x1b[nrow(x1b):1,])
s1$Line <- 1
s1 <- rbind(s1,s1[1,])
write.csv(s1,"stratum1.csv")
x11()
plot(s1[,1:2])
polygon(s1[,1:2],col=1)


points(x1b,col=2)
#### Stratum 2, grense i nord til 7 15', , i vest 80 30' (80.50)
x2 <- x1000[1264:1676,] ##Grense
x2b <- x20[1615:1937,]
s2 <- rbind(x2,x2b[nrow(x2b):1,])
s2$Line <- 2
s2 <- rbind(s2,s2[1,])
write.csv(s2,"stratum2.csv")
x11()
plot(s2)
polygon(s2,col=1)

points(x2,col=3)
#### Stratum 3, grense i vest 80 30' (80.50), grense i nord 6 20' (6.33)
x3 <- x1000[964:1263,] ##Grense
x3b <- x20[1324:1614,]
s3 <- rbind(x3,x3b[nrow(x3b):1,])
s3$Line <- 3
s3 <- rbind(s3,s3[1,])
write.csv(s3,"stratum3.csv")
x11()
plot(s3)
polygon(s3,col=1)
head(x3)
points(x3,col=4)

#### Stratum 4, grense i sør 6 20' (6.33), grense i nord 8 20' (8.33)
x4 <- x1000[544:963,] ##Grense
x4b <- x20[733:1323,]
s4 <- rbind(x4,x4b[nrow(x4b):1,])
s4$Line <- 4
s4 <- rbind(s4,s4[1,])
write.csv(s4,"stratum4.csv")
x11()
plot(s4)
polygon(s4,col=1)
head(x4)
points(x4,col=2)

#### Stratum 5, grense i sør 8 20' (8.33), grense i nord 9 30' (9.50)
x5 <- x1000[292:543,] ##Grense
x5b <- x20[390:732,]
s5 <- rbind(x5,x5b[nrow(x5b):1,])
s5$Line <- 5
s5 <- rbind(s5,s5[1,])
write.csv(s5,"stratum5.csv")
x11()
plot(s5)
polygon(s5,col=1)
head(x5)
points(x5,col=3)

#### Stratum 6, grense i sør 9 30' (9.50), grense i nord mot EZZ 10 30' ved 1000m (10.50)
x6 <- x1000[107:291,] ##Grense
x6b <- x20[208:389,]
s6 <- rbind(x6,x6b[nrow(x6b):1,])
s6$Line <- 6
s6 <- rbind(s6,s6[1,])
write.csv(s6,"stratum6.csv")
plot(s6)
polygon(s6)
head(x6)
points(x5,col=5)

