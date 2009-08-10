## ssm.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:34:41 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Mon Aug 10 15:32:29 2009
## Update Count    : 92
## Status          : Unknown, Use with caution!
###############################################################################

# source("f:/Research/sspir/sspirdevel/sspir/R/ssm.R")

getFit <- function(ssm) ssm$ss

m0.ssm <- function(ssm) ssm$ss$m0
C0.ssm <- function(ssm) ssm$ss$C0
Fmat.ssm <- function(ssm) ssm$ss$Fmat
Gmat.ssm <- function(ssm) ssm$ss$Gmat
Vmat.ssm <- function(ssm) ssm$ss$Vmat
Wmat.ssm <- function(ssm) ssm$ss$Wmat
phi.ssm <- function(ssm) ssm$ss$phi

"m0<-.ssm" <- function(ssm,value) {ssm$ss$m0<-value;return(ssm)}
"C0<-.ssm" <- function(ssm,value) {ssm$ss$C0<-value;return(ssm)}
"Fmat<-.ssm" <- function(ssm,value) {ssm$ss$Fmat<-value;return(ssm)}
"Gmat<-.ssm" <- function(ssm,value) {ssm$ss$Gmat<-value;return(ssm)}
"Vmat<-.ssm" <- function(ssm,value) {ssm$ss$Vmat<-value;return(ssm)}
"Wmat<-.ssm" <- function(ssm,value) {ssm$ss$Wmat<-value;return(ssm)}
"phi<-.ssm" <- function(ssm,value) {ssm$ss$phi<-value;return(ssm)}

"ssm" <-
  function(formula,
           family=gaussian,
	   data = list(),
           subset=NULL,
           fit=TRUE,
           phi=NULL,
           m0=NULL,
           C0=NULL,
           Fmat=NULL,
           Gmat=NULL,
           Vmat=NULL,
           Wmat=NULL
           )
{

  ## ###########################
  ## Helper functions

  noTvar <- function(e) {
    ## remove tvar() from formula
    if (is.call(e))
      if (e[[1]]==as.name("tvar")) 
        return(e[[2]])
      else {
        for (i in 2:length(e)) 
          e[[i]] <- noTvar(e[[i]])
        return(e)
      }
    else
      e
  }

  removeSpace <-  function(string) gsub("[[:space:]]", "", string)

  makeSS <-
    function(x,y,tvar,tvarNames,...) {
      
      dt = rep(1,nrow(x))
      p <- ncol(x)
      
      ntvar <- length(unique(tvar[tvar!=0]))
      phi <- c(1,rep(1,ntvar))
      names(phi) <- c("epsilon",tvarNames)
      
      res <- SS(
                y=ts(matrix(y,ncol=1),start=start(y),end=end(y),frequency=frequency(y)),
                x=list(x=x,dt=dt,tvar=tvar),
                phi = phi,
                Fmat = function(tt,x,phi) { x$x[tt,] },
                Gmat = function(tt,x,phi) { diag(p) },
                Vmat = function(tt,x,phi) { matrix(phi["epsilon"],1,1) },
                Wmat = function(tt,x,phi) {
                  W <- matrix(0,p,p)
                  diag(W)[tvar!=0] <- phi[-1][tvar]
                  return(W)
                },
                C0 = diag(p)*1e6,
                m0 = matrix(0,1,p)
                )

      return(res)
    }


  
  ## ###########################
  ## Handle the call (close to lm/glm)
  cl <- match.call()
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("`family' not recognized")
  }
  
  if (missing(data)) 
    data <- environment(formula)

  mf <- match.call(expand.dots = FALSE)

  ## ####################################################################
  ## Creates formula.notvar where tvar is removed from
  ## the formula

  ## special cases to be treated
  ## polytime
  ## season
  ## sumseason
  
  formula.notvar <- noTvar(formula)
  
  ## ####################################################################
  ## evaluate variables in the call
  mf$formula <- formula.notvar
  mf$family <- mf$start <- NULL
#  mf$... <- NULL
  mf$drop.unused.levels <- TRUE
  mf$na.action <- "na.pass"
  mf$fit  <- NULL
  mf$phi  <- NULL
  mf$m0  <- NULL
  mf$C0  <- NULL
  mf$Fmat  <- NULL
  mf$Gmat  <- NULL
  mf$Wmat  <- NULL
  mf$Vmat  <- NULL
  mf[[1]] <- as.name("model.frame")
  
#  mf <- eval(mf, parent.frame())
#  y <- eval(mf[[2]][[2]],.GlobalEnv)
 y <- eval(mf[[2]][[2]],data)  
  time <- time(ts(start=start(y),end=end(y),frequency=frequency(y)))
  assign("time",time,env=.GlobalEnv)
  mf <- eval(mf, .GlobalEnv)
  mt <- attr(mf, "terms")
    
#  na.act <- attr(mf, "na.action")
  xvars <- as.character(attr(mt, "variables"))[-1]
  if ((yvar <- attr(mt, "response")) > 0) 
    xvars <- xvars[-yvar]
  xlev <- if (length(xvars) > 0) {
    xlev <- lapply(mf[xvars], levels)
    xlev[!sapply(xlev, is.null)]
  }
  
  y <- model.response(mf, "numeric")
  w <- NULL
#  y <- ts(y,start=start,end=end,frequency=frequency)
#  w <- model.weights(mf)
#  offset <- model.offset(mf)

#  if (!is.null(weights) && any(weights < 0)) 
#    stop("Negative wts not allowed")
  
#  if (!is.null(offset) && length(offset) != NROW(y)) 
#    stop("Number of offsets is ", length(offset), ", should equal ", 
#         NROW(y), " (number of observations)")


  ## ##############################
  ## create model.matrix
  x <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf)
  else
    matrix(, NROW(y), 0)

  ## ########################################
  ## Figure out which terms are time-varying

  term.labels <- unname(sapply(attr(terms(formula.notvar), "term.labels"),
                               removeSpace))

  pAssign <- attr(x, "assign")
  if (attr(terms(formula),"intercept")) {
    pAssign <- pAssign + 1
    pAssign[1] <- 1
    term.labels <- c("(Intercept)",term.labels)
  }
    
  tvar.terms <- terms( formula, specials = "tvar",keep.order=TRUE )
  idx <- attr(tvar.terms,"specials")$tvar
  if (attr(tvar.terms,"response")) idx <- idx - 1

  if (attr(tvar.terms,"intercept")) 
    if (!any(attr(tvar.terms,"term.labels")=="tvar(1)"))
      idx <- idx + 1
  
  tvar <- attr(terms(formula.notvar,keep.order=TRUE),"term.labels")
  if (attr(tvar.terms,"intercept")) tvar <- c("(Intercept)",tvar)

  tvar <- tvar[idx]
  tvar <- unname( sapply( tvar, removeSpace) )

  tvarAssign <- match(pAssign, sort(match(tvar, term.labels)))
  tvarAssign[is.na(tvarAssign)] <- 0

  ## tvar now holds the names of the time-varying columns in the
  ## design matrix
  ## tvarAssign holds the indices (from 1 to No. of tvar() terms) of
  ## the time-varying columns so that 
  ## entries in the same tvar-statement has the same index. 

  ## ###############################
  ## Handle polytime
  ## ###############################
  if (any( unlist(lapply(term.labels,substr,1,8))=="polytime")) {
#    cat("polytime present\n")
    ## determine degree
    d <- sum(unlist(lapply(colnames(x),substr,1,8))=="polytime")

    nm <- unlist(lapply(colnames(x),substr,1,8))
    polytime.idx <-  (1:ncol(x))[nm=="polytime"]

    Pmat <- function(deg) {
      res <- matrix(0,deg+1,deg+1)
      diag(res) <- 1
      for (i in 1:deg) {
        res[1:i,i+1] <- 1/factorial(i:1)
      }
      res
    }

#    print(Pmat(d-1))
#    print(polytime.idx)

    ## if tvar(polytime, adjust tvar and tvarAssign
      if (any( unlist(lapply(tvar,substr,1,8))=="polytime")) {

        polnum <- pmatch("polytime",tvar)
#print(polnum)
#        print("---NU\n")
#        print(colnames(x)[polytime.idx])
#        print(tvar)
#        print("---NU\n")
#        print(tvarAssign)
        nytvarAssign <- tvarAssign
#        print(d)
        nytvarAssign[tvarAssign>polnum] <- tvarAssign[tvarAssign>polnum] + d-1
        nytvarAssign[tvarAssign==polnum] <- polnum:(polnum+d-1)
#        print(nytvarAssign)
        tvarAssign <- nytvarAssign

        if (polnum>1)
          nytvar <- tvar[1:(polnum-1)]
        else
          nytvar <- c()
        if (polnum<length(tvar))
          nytvar <- c(nytvar,colnames(x)[polytime.idx],tvar[ (polnum+1):length(tvar)])
#        print(nytvar)
        tvar <- nytvar
#        print("---NU\n")
#        print(tvar)
#        print(tvarAssign)
        
#        print("---NU\n")
      }
  }
    
  ## ###############################
  ## Handle sumseason
  ## ###############################
  if (any( unlist(lapply(term.labels,substr,1,9))=="sumseason")) {
    ## Count the number of terms including 'sumseason' in the column
    ## names of the design matrix. This is one less than the frequency
    period <- 1 + sum(unlist(lapply(colnames(x),substr,1,9))=="sumseason") 

    Cmat <- function(ncol) rbind( -1, cbind(diag(ncol-1),0) )

    nm <- unlist(lapply(colnames(x),substr,1,9))
    sumseason.idx <-  (1:ncol(x))[nm=="sumseason"]
    
  }
  
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = numeric(0), residuals = y,
              fitted.values = 0 * y, rank = 0, df.residual =
              length(y)) 
#    if (!is.null(offset)) 
#      z$fitted.values <- offset
  }
  else {
    ## if binomial, lhs is given as "cbind(y,n-y)"
    if (family$family=="binomial") {
      ntotal<- y[1,]+y[2,]
      y <- y[1,]
    }
    else {
      ntotal <- NA
    }

    ## ##############################
    ## Create State Space model
    z <- list()
    z$ss <- makeSS(x, y, w, tvar=tvarAssign,
                   tvarNames=tvar)

    if (any( unlist(lapply(term.labels,substr,1,9))=="sumseason")) {
      ## Adjust G-matrix if sumseason present
      oldG <- z$ss$Gmat
      newG <- function(tt,x,phi) {
        Gmat <- oldG(tt,x,phi)
        Gmat[sumseason.idx,sumseason.idx] <- Cmat(period-1)
        return(Gmat)
      }
      z$ss$Gmat <- newG

      ## adjust W-matrix
      if (any( unlist(lapply(tvar,substr,1,9))=="sumseason")) {
        oldW <- z$ss$Wmat
        newW <- function(tt,x,phi) {
          Wmat <- oldW(tt,x,phi)
          diag(Wmat)[sumseason.idx[-1]] <- 0
          return(Wmat)
        }
        z$ss$Wmat <- newW
      }
      
    }

    ## Adjust polytime
    if (any( unlist(lapply(term.labels,substr,1,8))=="polytime")) {
      oldG2 <- z$ss$Gmat
      newG2 <- function(tt,x,phi) {
        Gmat <- oldG2(tt,x,phi)
        Gmat[polytime.idx,polytime.idx] <- Pmat(d-1)
        return(Gmat)
      }
      z$ss$Gmat <- newG2
    }

  }
    class(z) <- c("ssm","lm")
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- xlev
    z$call <- cl
    z$terms <- mt
    z$data  <- data
    z$ss$family <- family
    z$ss$ntotal <- ntotal
  if (!is.null(phi)) z$ss$phi[1:length(phi)] <- phi
  if (!is.null(m0)) z$ss$m0 <- m0
  if (!is.null(C0)) z$ss$C0 <- C0
  if (!is.null(Fmat)) z$ss$Fmat <- Fmat
  if (!is.null(Gmat)) z$ss$Gmat <- Gmat
  if (!is.null(Vmat)) z$ss$Vmat <- Vmat
  if (!is.null(Wmat)) z$ss$Wmat <- Wmat
  if (fit) z$ss <- kfs(z)
  z
  
}

"print.ssm" <-
function(x,...) {
  print(x$ss$family)

  cat("Approximating SSM below:\n")
  print(x$ss)
}
