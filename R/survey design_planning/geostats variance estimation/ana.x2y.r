"ana.x2y" <-
  function (db,flag.plot=F) 
  {
    #===============================================================================
    # GAUSSIAN ANAMORPHOSIS
    #
    # UE program Fisboat, DG-Fish, STREP n° 502572
    # author: M.Woillez, ENSMP
    # date: 28 mars 2006 
    # Adapted by N. BEZ septembre 2013
    #
    # function: f.anamorphosis
    #           Compute the gaussian anamorphosis of the variable "z" of 
    #           the database "db". Create a gaussian variable named "gauss". 
    #===============================================================================
    
    ########## sub-routine fdr()
    
    fdr <- function (x,w,flag.plot=F)
    {
      # To perform a suitable cumulative distribution function
      # For every real number x, the cdf is given by:
      # F(x) = P[X<=x] or sum(P[X=x]) for X<=x
      # Here P[X=x] equals the area of influence of the data x divided by the 
      # total area of influence (e.g. the weights w).
      # Then P[X<=xk] is modified as pk=(pk + pk-1)/2
      # to avoid an indetermination for the gaussian quantile function when p=1 
      
      # sort the variable x and its weigths w 
      xx <- sort(x,index.return=T)
      x <- xx$x
      w <- w[xx$ix]
      
      # set the probability P[X<=x] using w 
      p<-numeric(0)
      for(i in 1:length(x))
        p[i]<-sum(w[1:i])
      
      # modify the probability as pk=(pk + pk-1)/2 to avoid an indetermination for the gaussian quantile function when p=1
      p<-((c(0,p[-length(p)])+p)/2)
      
      # reduce redondancy in the data x  
      # keep the associated cumulative probability
      # and count the number of equal values
      a<-1
      j<-1 
      xc<-x[1]
      b<-numeric(0)
      pc<-numeric(0)
      for(i in 1:(length(x)-1)){
        if(x[i+1]!=x[i]){
          xc<-c(xc,x[i+1])
          pc<-c(pc,p[i])
          j<-i+1
          b<-c(b,a)
          a<-1
        }
        else a<-a+1  
      }
      b<-c(b,a)
      pc<-c(pc,p[j])
      
      # plot the cumulative distribution
      if(flag.plot==T){
        x11()
        wind(2,2)
        plot(xc,pc,pch=".",col=1,xlab="Z",ylab="F(Z)",ylim=c(0,1))
        lines(xc,pc)
      }
      # output
      res<-list("z"=x,"class"=xc,"nb"=b,"prob"=p,"pcum"=pc,"rank"=xx$ix)
      res
    }
    ############## Fin de la sub-routine
    
    z <- db[,"z1"]
    w <- db[,"w"]
    w <- w/sum(w)
    
    # make the cumulative distribution and store the wanted results in new objects
    temp <- fdr(z,w)
    pcum <- temp$pcum
    nb <- temp$nb
    rank <- temp$rank
    
    # Generate the gaussian variable with the function "qnorm" which gives the quantile function
    y<- qnorm(pcum)
    
    # plot the cumulative distribution
    if(flag.plot==T){
      plot(y,pcum,pch=".",col=1,xlab="Y+",ylab="G(Y+)",ylim=c(0,1))
      lines(y,pcum)
    }
    
    # expand the gaussian variable according to the previous redondancy
    flag <- 0
    for(j in 1:length(y)){
      if(flag==0){
        bid <- rep(y[j],nb[j])
        flag <- 1
      }
      else bid <- c(bid, rep(y[j],nb[j]))
    }
    
    # rank the gaussian variable according to the original data x
    Y <- rep(0,length(bid))
    for(k in 1:length(bid)){
      Y[rank[k]] <- bid[k]
    }
    
    # set the output in a correct format for the database db
    # the new variable is named "gauss" and is active with the locator "z" 
    # in the database db
    res <- rep(NA,dim(db@items)[1])        
    #res[data$valid.index]<-Y
    db <- db.add(db, "gauss" = Y)
    db <- db.locate(db, "gauss", "z",start=1)
    db
  }
