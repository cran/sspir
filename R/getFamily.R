## getFamily.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:29:59 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Mon Aug 10 15:31:46 2009
## Update Count    : 7
## Status          : Unknown, Use with caution!
###############################################################################

"getFamily" <-
function(family) {
  ## family object contains link and family
  if (family$family=="poisson") {
    ## general poisson
    initV <- function(y) { rep(1,length(y)) }
    
    initY <- function(y) { res <- log(y);res[y==0] <- -1;res}
    simY   <- function(n,lambda,ntotal,...) rpois(n,family$linkinv(lambda))
                         
    ## specific poisson-link combination
    if (family$link=="log") {
      vtilde <- function(lambda,ntotal,...)  exp(-lambda)	
      ytilde <- function(lambda,y,vtilde,...) lambda + vtilde*y - 1
    }
    else stop("Extended Kalman not implemented for ",family$family,family$link," combination")
    
  }
  else if (family$family=="binomial") {
    ## general binomial
    initV <- function(y) { rep.int(1,length(y))
    }
    initY <- function(y) expit(y)
    simY   <- function(n,lambda,ntotal,...) rbinom(n,ntotal,family$linkinv(lambda))

    ## specific binomial-link combination
    if (family$link=="logit") {
      vtilde <- function(lambda,ntotal,...) {
        res <- ((1+exp(lambda))^2)/(ntotal*exp(lambda))
        res[ntotal==0] <- 1/4
        return(res)
      }
      ytilde <- function(lambda,y,vtilde,...) lambda + vtilde*y - (1+exp(lambda))
    }
    else stop("Extended Kalman not implemented for ",family$family,family$link," combination")

    
  }
  else stop("Extended Kalman not implemented for ",family$family,family$link," combination")

  list(initV=initV,initY=initY,vtilde=vtilde,ytilde=ytilde,simY=simY)
}

"logit" <-
function(mu) log(mu/(1 - mu))

"expit" <-
function(x) exp(x)/(1+exp(x))
