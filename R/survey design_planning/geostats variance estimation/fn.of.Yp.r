"fn.of.Yp"<-
function (yc, n, norm = T)
{
#===============================================================================
# Coefficient du développement de la fonction f(Y)=Yp en polynomes d'Hermite
# fn[i+1] correspond à fn_i
# fn_0 = fn[1] est fourni même s'il doit être supprimé ensuite dans le lien entre 
# cov(Yp) et cov(Y)

H <- hermite(y=yc, n=n, norm = T)
fn <- dnorm(yc) + yc*pnorm(yc)	### fn_0
fn[2] <- pnorm(yc)-1			### fn_1
for(i in 3:n){
	fn[i] <- (yc*dnorm(yc)*H[i-1]/sqrt(i-1)) + (H[i-2]*dnorm(yc)/sqrt((i-1)*(i-2)))
}
fn
}
    




