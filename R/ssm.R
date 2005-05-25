## ssm.R --- 
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:34:41 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Sun Jan 23 16:03:53 2005
## Update Count    : 8
## Status          : Unknown, Use with caution!
###############################################################################

"ssm" <-
  function(formula,
           family=gaussian,
	   data = list(),
           subset=NULL,
           time=NULL
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

  ## makeSS.R --- 

"makeSS" <-
function(x,y,tvar,tvarNames,...) {
#  printline()
#  cat("makeSS\n")
#  cat("x:\n");print(x)
#  cat("y:\n");print(y)
#  cat("w:\n");print(w)
  
#  cat("tvar:\n");print(tvar)
  
#  if (is.null(w)) weights <- FALSE
#  else warning("Cannot handle weights, yet.\n")

  dt = rep(1,nrow(x))
  p <- ncol(x)

  ntvar <- length(unique(tvar[tvar!=0]))
  phi <- c(1,rep(1,ntvar))
  names(phi) <- c("epsilon",tvarNames)
  
  res <- SS(
            y=matrix(y,nrow=1),
            x=list(x=x,dt=dt,tvar=tvar),
            phi = phi,
            Fmat = function(tt,x,phi) { x$x[tt,] },
            Gmat = function(tt,x,phi) { diag(p) },
            Vmat = function(tt,x,phi) { matrix(phi["epsilon"],1,1) },
            Wmat = function(tt,x,phi) {
              W <- matrix(0,p,p)
              diag(W)[tvar!=0] <- phi[-1][tvar]
              return(W)
              ## Her skal ganges med dt
            },
##            C0 = matrix(var(y)*1e4,p,p),
            C0 = diag(p)*1e6,
            m0 = matrix(0,p,1)
            )

  #list(x,y,w,offset,tvar)
  return(res)
}


  
  ## ###########################
  
  
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
  mf$time  <- NULL
  mf[[1]] <- as.name("model.frame")
  
  mf <- eval(mf, parent.frame())
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
#  w <- model.weights(mf)
#  offset <- model.offset(mf)

#  if (!is.null(weights) && any(weights < 0)) 
#    stop("Negative wts not allowed")
  
#  if (!is.null(offset) && length(offset) != NROW(y)) 
#    stop("Number of offsets is ", length(offset), ", should equal ", 
#         NROW(y), " (number of observations)")

  if (missing(time)) {
    time <- 1:length(y)
    warning("time set to 1:length(y)\n")
  }

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
      ntotal<- y[,1]+y[,2]
      y <- y[,1]
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

#    if (length(pf$terms.notvar)>0&&pf$terms.notvar=="") {
#        ## model: y~tvar(). Trend model.
#        z$ss$Wmat <- function(tt,x,phi) { matrix(c(phi[1],0,0,phi[2]),2,2) }
#        z$ss$Vmat <- function(tt,x,phi) { matrix(phi[3],1,1) }
#        z$ss$Gmat <- function(tt,x,phi) { matrix(c(1,0,1,1),2,2) }
#        z$ss$Fmat <- function(tt,x,phi) { matrix(c(1,0),2,1) }
#        z$ss$phi  <- c(1,1)
#        z$ss$p    <- 2
#        z$ss$m0   <- matrix( 0, 2,1)
#        z$ss$C0   <- diag(2)*1e6
#      }

    }
    class(z) <- c("ssm","lm")
#    if (!is.null(na.act)) 
#        z$na.action <- na.act
#    z$offset <- offset
    z$contrasts <- attr(x, "contrasts")
    z$xlevels <- xlev
    z$call <- cl
    z$terms <- mt
    z$data  <- data
    z$ss$family <- family
    z$ss$ntotal <- ntotal
    z
}

