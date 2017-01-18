"ana.model" <-
function(db,name.Z,name.Yp)
{
par(mfrow=c(2,2))
Yp <- db[,name.Yp]
Zp <- db[,name.Z]
hist(Zp,nclass=100,main="Histogram of Z",xlab="",ylab="",col=8)
hist(Yp,nclass=50,main="Histogram of Y+",xlab="",ylab="",col=8)
ycut<-min(Yp,na.rm=T)
y<-Yp[Zp>0]
z<-Zp[Zp>0]
y<-y[!is.na(y) & !is.na(z)]
z<-z[!is.na(y) & !is.na(z)]
yy<-y[rev(order(y))]
zz<-z[rev(order(y))]
out<-lm(zz[1:2]~yy[1:2])
yout <- 5
y<-c(y,yout)
z<-c(z,out$coeff[2]*yout+out$coeff[1])
plot(y,z,xlab="Y+>ycut",ylab="Z>0",main="Gaussian anamorphosis",xlim=c(ycut,max(yy)),ylim=c(0,max(zz)))
f <- approxfun(y,z,rule=2)
curve(f(x),ycut,max(y),col=2,add=T)
qqnorm(Yp)
rm(Yp,Zp,y,z,yy,zz,out,yout)
par(mfrow=c(1,1))
return(f)
}