## SS.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:34:35 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Sun Jan 23 16:03:17 2005
## Update Count    : 4
## Status          : Unknown, Use with caution!
###############################################################################

SS <- function(
               y=NA,
               x=NA,
               Fmat = function(tt,x,phi) { NA },
               Gmat = function(tt,x,phi) { NA },
               Vmat = function(tt,x,phi) { NA },
               Wmat = function(tt,x,phi) { NA },
               m0 = 0,
               C0 = NA,
               phi = NA
               ) {

  ## Arguments
  ## y: data vector (optional, since it can be simulated by 'recursion')
  ## x: list of covariates (optional)
  ## Fmat(tt,x,phi): function to give design matrix at time tt depending on x and phi
  ## Gmat(tt,x,phi): function to give innovation matrix at time tt depending on x and phi
  ## Vmat(tt,x,phi): function to give obs. variance matrix at time tt depending on x and phi
  ## Wmat(tt,x,phi): function to give state. variance matrix at time tt depending on x and phi
  ## m0: Starting place for state vector
  ## C0: Variance of m0
  ## phi: parameter-vektor

  ss <- list(y=y, x=x, Fmat=Fmat, Gmat=Gmat, Vmat=Vmat, Wmat=Wmat, m0=m0, C0=C0, phi=phi)

  ## the function sets up an SS object with attributes
  ## n: length of y (if y is given)
  ## d: dimension of y (if y is given). 
  ## p: dimension of theta

  if (length(y)==1) {
    if (is.na(y))
      n <- NA
    else n <- 1
    d <- dim(ss$Fmat(1,ss$x,ss$phi))[2]
  }
  else
    {
      n <- dim(y)[2]
      d <- dim(y)[1]
      ##    if (length(dim(d))==0) d <- 1
    }

  if ( length(ss$Gmat(1,ss$x,ss$phi))>1) {
    p <- dim(ss$Gmat(1,ss$x,ss$phi))[1]
  }
  else
    {
    if ( length(ss$Gmat(1,ss$x,ss$phi)) == 1) p <- 1
    else
      p <- NA
  }
  
  ss <- c(ss, list(n=n, d=d, p=p) )
  
  ## ytilde: used in extended 
  ## iteration: Number of iterations in extended
  ## m: E[theta|y's]. Without m0.
  ## C: Var[theta|y's]. Without C0.
  ## mu: F^t * theta. Signal
  ## likelihood: loglikelihood( phi )

  ytilde     <- NA
  iteration  <- 0
  m          <- NA
  C          <- NA
  mu         <- NA
  likelihood <- NA

  ss <- c(ss, list(ytilde=ytilde, iteration=iteration,
                   m=m,C=C,mu=mu,likelihood=likelihood) )
  
  class(ss) <- "SS"
  ss
}

## print.SS.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:33:25 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Fri Jan 21 18:35:31 2005
## Update Count    : 2
## Status          : Unknown, Use with caution!
###############################################################################

"print.SS" <-
function(x,...) {

  printline()
  cat("The state space model is given by\n")
  cat("\n")
  cat("Y_t     = F_t^T %*% theta_t     + v_t, v_t ~ N(0,V_t)\n")
  cat("theta_t = G_t   %*% theta_{t-1} + w_t, w_t ~ N(0,W_t)\n")
  cat("\n")
  cat("for t=1,...,",x$n,"\n",sep="")
  printline()

  cat("Dimensions of the involved terms:\n")

  cat("dim(Y_t) =", x$d, "x 1   (dx1)\n")
  cat("dim(F_t) =", x$p, "x",x$d,"  (pxd)\n" )
  cat("dim(G_t) =", x$p, "x",x$p,"  (pxp)\n" )
  cat("dim(V_t) =", x$d, "x", x$d,"  (dxd)\n" )
  cat("dim(W_t) =", x$p, "x",x$p,"  (pxp)\n" )
  cat("dim(C_0) =", x$p, "x",x$p,"  (pxp)\n" )
  cat("dim(m_0) =", x$p, "x 1   (px1)\n" )
  
  printline()

  cat("The terms:\n")

  cat("Y=\n")
  if (length(x$y)<10)
    print(x$y)
  else {
    cat("First 10 cols\n")
    print(x$y[,1:10])
  }
  
  cat("phi= (the parametervector)\n")
  print(x$phi)

  cat("\nx= (the list of covariates)\n")
  if (length(x$x$x)<10)
    print(x$x)
  else {
    cat("First 10 rows\n")
    print(x$x$x[1:10,])
  }

  cat("\nFormula for creating F_t:\n")
  print(x$Fmat)
  cat("\nF_1^T = \n")
  print(t(x$Fmat(1,x$x,x$phi)))

  cat("\nFormula for creating V_t:\n")
  print(x$Vmat)
  cat("\nV_1 = \n")
  print(x$Vmat(1,x$x,x$phi))

  cat("\nFormula for creating G_t:\n")
  print(x$Gmat)
  cat("\nG_1 = \n")
  print(x$Gmat(1,x$x,x$phi))
  
  cat("\nFormula for creating W_t:\n")
  print(x$Wmat)
  cat("\nW_1 = \n")
  print(x$Wmat(1,x$x,x$phi))

  cat("\nm0 = \n")
  print(x$m0)

  cat("\nC0 = \n")
  print(x$C0)
  
}

## plot.SS.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:32:43 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Fri Jan 21 12:32:44 2005
## Update Count    : 1
## Status          : Unknown, Use with caution!
###############################################################################

"plot.SS" <-
function(x,...) {
  ss <- x
  
  xs <- 1:ss$n

  plot(xs, ss$y)
  
  lines(xs, ss$mu, lty=2)

  if (ss$p == 1) {
    lines(xs, ss$mu - 2*sqrt(unlist(ss$C)), lty=3)
    lines(xs, ss$mu + 2*sqrt(unlist(ss$C)), lty=3)
  }

  if (length(ss$truetheta)>0 && ss$p == 1)
    {
      points(xs, ss$truetheta,col="dark red",pch=3)
      cat("True theta are marked with dark red crosses\n")
    }
  

  invisible()
}

## print.Smoothed.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:33:17 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Fri Jan 21 12:33:19 2005
## Update Count    : 1
## Status          : Unknown, Use with caution!
###############################################################################

"print.Smoothed" <-
function(x,...) {
  print.SS(x)

  printline()

  cat("After smoothing...\n")

  cat("Mu = (the signal, F_t^T m_t) \n")
  print(x$mu)
  printline()
  
  cat("log(Likelihood) =", x$likelihood,"\n")

  cat("(Note that m0 and C0 has been replaced by E(m0|Y) and E(C0|Y) )\n")
  printline()
  
  cat("m_1, m_2, m_3 =\n")
  print(x$m[,1:3])

  cat("C_1, C_2, C_3 =\n")
  print(x$C[[1]])
  print(x$C[[2]])
  print(x$C[[3]])
  invisible()
}

