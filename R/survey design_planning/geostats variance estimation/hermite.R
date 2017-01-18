"hermite" <-
  function (y, n, norm = T)
  {
    #===============================================================================
    # Polynomes d'Hermite
    # H[i+1] correspond à H_i
    # H_0=H[1]=1 etc
    
    if(norm == T){
      H <- rep(NA,n)
      H[1] <- 1
      H[2] <- -y
      H[3] <- (1/sqrt(2))*(y^2-1)
      if(n >= 4){
        for(i in 2:(n-2)){
          H[i+2] <- -(1/sqrt(i+1))*y*H[i+1]-sqrt(i/(i+1))*H[i]
        }
      }
    }
    else{
      H <- rep(0,n)
      H[1] <- 1
      H[2] <- -y
      H[3] <- (y^2-1)
      if(n >= 4){
        for(i in 2:(n-2)){
          H[i+2] <- (-y*H[i+1]) -(i* H[i])
        }
      }
    }
    H
  }
