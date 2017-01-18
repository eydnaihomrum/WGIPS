"fn.of.Yp"<-
function (yc, n, norm = T)
{
#===============================================================================
# Coefficient du d�veloppement de la fonction f(Y)=Yp en polynomes d'Hermite
# fn[i+1] correspond � fn_i
# fn_0 = fn[1] est fourni m�me s'il doit �tre supprim� ensuite dans le lien entre 
# cov(Yp) et cov(Y)

H <- hermite(y=yc, n=n, norm = T)
fn <- dnorm(yc) + yc*pnorm(yc)	### fn_0
fn[2] <- pnorm(yc)-1			### fn_1
for(i in 3:n){
	fn[i] <- (yc*dnorm(yc)*H[i-1]/sqrt(i-1)) + (H[i-2]*dnorm(yc)/sqrt((i-1)*(i-2)))
}
fn
}
    




