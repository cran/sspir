## polytime.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:33:01 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Tue Apr 11 11:41:41 2006
## Update Count    : 3
## Status          : Unknown, Use with caution!
###############################################################################

polytime <- function(time,degree=1) {
  n <- length(time)
  p <- degree
  res <- matrix( c(1,rep(0,p)), n,p+1,byrow=TRUE)
  colnames(res) <- paste("time",0:p,sep="")
  res	
}

## polytrig.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:33:09 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Fri Jan 21 12:33:10 2005
## Update Count    : 1
## Status          : Unknown, Use with caution!
###############################################################################

"polytrig" <-
function(time,period=365.25,degree=1) {
    k <- 2*pi/period
    i <- 1
    res <- cbind(cos(time* i*k),sin(time* i*k))
    lab <- c("c1","s1")
    if (degree>1) {
        for (i in 2:degree) {
            res <- cbind(res,cos(time* i*k),sin(time* i*k))
            lab <- c(lab, paste("c",i,sep=""),paste("s",i,sep=""))
        }
    }
    colnames(res) <- lab
    res
}

## season.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:34:07 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Fri Jan 21 12:34:09 2005
## Update Count    : 1
## Status          : Unknown, Use with caution!
###############################################################################

"season" <-
function(time,period=12) {
#    warning("saeson not implemented\n")

    position <-  time %% period 
    position[position==0] <- period

    value <- factor(position,labels=paste("p",1:period,sep=""))
    
#    value <- matrix(0,nrow=length(time),ncol=period)
#    for (i in 1:nrow(value)) value[i,position[i]] <- 1

#    colnames(value) <- paste("p",1:period,sep="")
    
    return(value)
}

sumseason <- function(time,period=12) {
  ## As in StructTS, the sum of contributions of the last 'period'
  ## terms are zero (or white noise if surrounded by tvar).

  ## How to handle time gaps? (force user to insert missing observations)
#  if (!all(diff(time)==diff(time)[1]))
#    warning("Not handling time gaps")
  
  
  ## G-matrix contribution:
  # Cmat <- function(ncol) rbind( -1, cbind(diag(ncol-1),0) )

  ## F-matrix contribution
  Fmat <- function(ncol) c(1,cbind(rep(0,ncol-2)))
  
  matrix(Fmat(period),nrow=length(time),ncol=period-1,byrow=TRUE)
}
