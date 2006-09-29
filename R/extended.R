## extended.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:30:11 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Thu Apr 20 13:00:38 2006
## Update Count    : 7
## Status          : Unknown, Use with caution!
###############################################################################

"extended" <-
function(ss,maxiter=50,epsilon=1e-6,debug=FALSE) {
  UseMethod("extended")
}

"extended.SS" <-
function(ss,maxiter=50,epsilon=1e-6,debug=FALSE)
{
  ## Run smoother(kfilter(ss)) iteratively, in each step adjusting the
  ## model from Non-gaussian to Gaussian. Return m, C, m0, C0, x

  ## ytilde and x$wtilde and x$vtilde are assumed initialized before the call to extended.

  converged <- FALSE
  iter <- 0
  origy <- ss$y
  origC0<- ss$C0
  origm0<- ss$m0

  fam <- getFamily(ss$family)
  ss$ytilde   <- fam$initY(ss$y)
  ss$x$vtilde <- fam$initV(ss$y)
  ss$Vmat <- function(tt,x,phi) { matrix(ss$x$vtilde[tt],1,1) }
  
  while (!converged & iter < maxiter) {

    oldv <- ss$x$vtilde
    oldw <- ss$x$wtilde
    oldy <- ss$ytilde
    ss$y <- ss$ytilde

    ss <- smoother(kfilter(ss))

    ss$C0 <- origC0
    ss$m0 <- origm0
    ss$x$vtilde <- fam$vtilde(ss$mu,ss$ntotal)
    ss$ytilde   <- fam$ytilde(ss$mu,origy,ss$x$vtilde)
    vconv <- abs( (ss$x$vtilde - oldv)/oldv )
 
    conv <- max( vconv,  (ss$ytilde-oldy)/oldy)

    
    iter <- iter + 1
    converged <- conv < epsilon
    if (debug) cat("(",iter,") log l(",ss$phi,") = ", ss$loglik,"(conv=",conv,")\n")
  }

  ss$iteration <- iter
  ss$y <- origy
  cat("\n")
  class(ss) <- c("Extended","SS")
  ss
}
