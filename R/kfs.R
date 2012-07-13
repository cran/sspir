## kfs.R ---
## Author          : Claus Dethlefsen
## Created On      : Fri Jan 21 12:31:36 2005
## Last Modified By: Claus Dethlefsen
## Last Modified On: Thu Nov 19 08:43:59 2009
## Update Count    : 42
## Status          : Unknown, Use with caution!
###############################################################################

"kfs" <-
function(ss,...) {
  UseMethod("kfs")
}

kfs.ssm <- function(ss,...) {
  if (ss$ss$family$family=="gaussian")
    return(smoother(kfilter(ss$ss)))
  else
    return(extended(ss$ss))
}

kfs.SS <- function(ss,...)
  smoother(kfilter(ss))

"recursion" <-
function(ss,n) {
  UseMethod("recursion")
}

"recursion.SS" <-
function(ss,n) {

  ss$n <- n

  p <- ss$p
  d <- ss$d

  m0 <- t(ss$m0)
  C0 <- ss$C0

  Gmat <- ss$Gmat
  Fmat <- ss$Fmat
  Vmat <- ss$Vmat
  Wmat <- ss$Wmat

  phi <- ss$phi
  x   <- ss$x

#  theta <- matrix(NA,p,n)
#  y     <- matrix(NA,d,n)
#  mu    <- matrix(NA,d,n)
  theta <- matrix(NA,n,p)
  y     <- matrix(NA,n,d)
  mu    <- matrix(NA,n,d)

if(class(Wmat)=="matrix"){
  epswrand <- mvrnorm(n,rep(0,p),Wmat)
} else{
  epswrand <- mvrnorm(n,rep(0,p),Wmat(1,x,phi))
}


if(class(Vmat)=="matrix"){
  epsvrand <- mvrnorm(n,rep(0,d),Vmat)
} else {
epsvrand <- mvrnorm(n,rep(0,d),Vmat(1,x,phi))
}

  theta[1,] <- Gmat(1,x,phi) %*% m0 + epswrand[1,]
  mu[1,]    <- t(Fmat(1,x,phi)) %*% theta[1,]

    y[1,]     <-  mu[1,] + epsvrand[1,]
#  else
#      y[,1] <- fam$simY(d,mu[,1],ntotal)

  for (tt in 2:n) {
    theta[tt,] <- Gmat(tt,x,phi) %*% theta[tt-1,] + epswrand[tt,]
    mu[tt,]    <- t(Fmat(tt,x,phi)) %*% theta[tt,]
      y[tt,]     <-  mu[tt,] + epsvrand[tt,]

  }
  ss$y <- y
  ss$truetheta <- theta
  ss$mu <- mu

  ss
}



"kfilter" <-
function(ss) {
  UseMethod("kfilter")
}

"filterstep" <-
function(y,Fmat,Gmat,Vt,Wt,mx,Cx)
  {
    ## do one step in the Kalman filter (see kfilter for explanations)
    ## The inverse filter should give the same
    ## result as the ordinary filter, but be faster if the state vector is
    ## smaller than the obs.vector
    ## the dimension of Q is dxd and the dimension of R is qxq
    a <- Gmat %*% t(mx)
    R <- Gmat %*% Cx %*% t(Gmat) + Wt
#    class(R) <- c("Hermitian","matrix") # p.d.

    f <- t(Fmat) %*% a
    Q <- t(Fmat) %*% R %*% Fmat + Vt

        e <- y - f
        A <- R %*% Fmat %*% mysolve(Q)
        m <- a + A%*%e
        C <- R - A%*%Q%*%t(A)

#    if (length(y)>1) loglikterm <-log(dmvnorm(as.numeric(y),as.numeric(f),Q))
#    else             loglikterm <- -0.5*( log(2*pi) + log(Q) + (y-f)^2/Q)

    #if (length(y)>1){
        loglikterm <- dmvnorm(as.numeric(y),as.numeric(f),Q,log=TRUE)
    #} else {
    #    loglikterm <- -0.5*( log(2*pi) + log(Q) + (y-f)^2/Q)
    #}

    assign(".Last.m", m, envir=.GlobalEnv)
  list(m=m,C=C,loglikterm=loglikterm)
  }



"kfilter.SS" <-
function(ss) {

  ## Observation eq:    y[,t] = Fmat^T * theta_t  + N_d(0,diag(V[,t]))
  ## Evolution   eq:  tht_t   = Gmat * theta_{t-1}+ N_q(0,diag(W[,t]))
  ## Initial       :  tht_0   ~ N(m0,C0)
  ## The Kalman filter yields m_t and C_t, with
  ##            (tht_t|y[,1:t]) ~ N(m_t, C_t)
  ##
  ## In the general case, the recursion goes:
  ##  a_t = G_t * m_{t-1}
  ##  R_t = G_t C_{t-1}*G_t^T + W_t
  ##  f_t = F_t^T * a_t
  ##  Q_t = F_t^T * R_t * F_t + V_t
  ##  e_t = y_t - f_t
  ##  A_t = R_t * F_t * Q_t^{-1}
  ##  m_t = a_t + A_t * e_t
  ##  C_t = R_t - A_t * Q_t * A_t^T
  ##
  ## in our case, F_t = Fmat, G_t = Gmat, V_t = diag(V[,t]),
  ##              W_t = diag(W[,t])

#  m <- matrix(NA,ss$p,ss$n)
  m <- matrix(NA,ss$n,ss$p)
  C <- vector("list",ss$n)
if(class(ss$Vmat)=="matrix" & class(ss$Wmat)=="matrix")
	{
  firststep <-
    filterstep(
               matrix(ss$y[1,]),
               ss$Fmat(1,ss$x,ss$phi),
               ss$Gmat(1,ss$x,ss$phi),
               ss$Vmat,
               ss$Wmat,
               ss$m0,
               ss$C0
               )
	} else {  if(class(ss$Vmat)=="matrix" & class(ss$Wmat)!="matrix")
	           {
              firststep <-
              filterstep(
               matrix(ss$y[1,]),
               ss$Fmat(1,ss$x,ss$phi),
               ss$Gmat(1,ss$x,ss$phi),
               ss$Vmat,
               ss$Wmat(1,ss$x,ss$phi),
               ss$m0,
               ss$C0
               )
	           } else {  if(class(ss$Vmat)!="matrix" & class(ss$Wmat)=="matrix")
	                     {
                         firststep <-
                            filterstep(
                                matrix(ss$y[1,]),
                                ss$Fmat(1,ss$x,ss$phi),
                                ss$Gmat(1,ss$x,ss$phi),
                                ss$Vmat(1,ss$x,ss$phi),
                                ss$Wmat,
                                ss$m0,
                                ss$C0
                                )
	                     } else {
                                firststep <-
                                  filterstep(
                                  matrix(ss$y[1,]),
                                  ss$Fmat(1,ss$x,ss$phi),
                                  ss$Gmat(1,ss$x,ss$phi),
                                  ss$Vmat(1,ss$x,ss$phi),
                                  ss$Wmat(1,ss$x,ss$phi),
                                  ss$m0,
                                  ss$C0
                                  )
	                             }
                    }
            }
  m[1,] <- firststep$m
  C[[1]]<- firststep$C
  loglik<- firststep$loglikterm

  ## run the recursion
  for (tt in 2:ss$n)
    {
#       cat(tt," ", sep="")
if(class(ss$Vmat)=="matrix" & class(ss$Wmat)=="matrix")
	{
      nextstep <-
        filterstep(
                   matrix(ss$y[tt,]),
                   ss$Fmat(tt,ss$x,ss$phi),
                   ss$Gmat(tt,ss$x,ss$phi),
                   ss$Vmat,
                   ss$Wmat,
                   matrix(m[tt-1,],nrow=1),
                   C[[tt-1]]
                   )
	}	else { if(class(ss$Vmat)=="matrix" & class(ss$Wmat)!="matrix")
	           {
              nextstep <-
                filterstep(
                   matrix(ss$y[tt,]),
                   ss$Fmat(tt,ss$x,ss$phi),
                   ss$Gmat(tt,ss$x,ss$phi),
                   ss$Vmat,
                   ss$Wmat(tt,ss$x,ss$phi),
                   matrix(m[tt-1,],nrow=1),
                   C[[tt-1]]
                   )
              } else { if(class(ss$Vmat)!="matrix" & class(ss$Wmat)=="matrix")
	                     {
                          nextstep <-
                              filterstep(
                                  matrix(ss$y[tt,]),
                                  ss$Fmat(tt,ss$x,ss$phi),
                                  ss$Gmat(tt,ss$x,ss$phi),
                                  ss$Vmat(tt,ss$x,ss$phi),
                                  ss$Wmat,
                                  matrix(m[tt-1,],nrow=1),
                                  C[[tt-1]]
                                  )
                        } else {
                                nextstep <-
                                  filterstep(
                                    matrix(ss$y[tt,]),
                                    ss$Fmat(tt,ss$x,ss$phi),
                                    ss$Gmat(tt,ss$x,ss$phi),
                                    ss$Vmat(tt,ss$x,ss$phi),
                                    ss$Wmat(tt,ss$x,ss$phi),
                                    matrix(m[tt-1,],nrow=1),
                                    C[[tt-1]]
                                    )
	                             }
	                 }
	       }

      m[tt,]  <- nextstep$m
      C[[tt]] <- nextstep$C
      loglik  <- loglik + nextstep$loglikterm
    }
  if (is.ts(ss$y))
    ss$m <- ts(m,start(ss$y),end=end(ss$y),frequency=frequency(ss$y))
  else
    ss$m <- m
  ss$C <- C
  ss$likelihood <- loglik
  ss$loglik <- loglik

  ss
}

kfilter.ssm <- function(ss,...) {
  if (ss$ss$family$family=="gaussian")
    return(kfilter(ss$ss))
  else {
    ## approximate
    res <- extended(ss$ss)
    origy <- res$y
    res$y <- res$ytilde
    kres <- kfilter(res)
    kres$ss$y <- origy
    return(kres)
  }
}

#####################################################################
## Functions using KFAS
## Modified by     : Anette Luther Christensen
## Modified on     : 15 March 2012
####################################################################
Fkfilter <- function(ss, tvar){
    ## Kalman filtering using KFAS
    ## Check of time-varying matrices:
    Zt <- sspir2kfas(ss$Fmat, dim=c(ss$d, ss$p, tvar[1]), x=ss$x, phi=ss$phi)
    Tt <- sspir2kfas(ss$Gmat, dim=c(ss$p, ss$p, tvar[2]), x=ss$x, phi=ss$phi)
    Ht <- array(sspir2kfas(ss$Vmat, dim=c(ss$d, ss$d, tvar[3]), x=ss$x, phi=ss$phi), dim=c(ss$d, ss$d, tvar[3]))
    Qt <- sspir2kfas(ss$Wmat, dim=c(ss$p, ss$p, tvar[4]), x=ss$x, phi=ss$phi)

    # Changing initial values to fit to KFAS
    a1 <- Tt[,,1] %*% matrix(ss$m0, ncol=1)
    P1 <- (Tt[,,1] %*% ss$C0 %*% t(Tt[,,1]) + Qt[,,1])

    m1 <- KFAS::SSModel(y=ss$y, Z=Zt, T=Tt, R=diag(ss$p), H=Ht, Q=Qt, a1=a1, P1=P1, distribution="Gaussian" )

    kfas.filt <- KFAS::KFS(m1, smoothing="state")

    # Backtracking filtered means, m, and covariances, C
    ss$mf <- t( (solve(Tt[,,1:tvar[2]]) %*% kfas.filt$a)[,-c(1)] )
    ss$Cf <- list()
    # Calculating log likelihood value
    loglik <- dmvnorm(x = ss$y[1,], mean = matrix( Zt[,,1], nrow = 1) %*% Tt[,,1] %*% matrix( ss$m0, ncol = 1), sigma = matrix( Zt[,,1], nrow = 1) %*% (Tt[,,1] %*% ss$C0 %*% t(Tt[,,1]) + Qt[,,1]) %*% matrix( Zt[,,1], ncol = 1) + Ht[,,1], log = TRUE)
    for(j in 1:ss$n){
        ss$Cf[[j]] <- solve(Tt[,,(tvar[2]!=1)*(j+2)+(tvar[2]==1)]) %*% (kfas.filt$P[,,j+1] - Qt[,,(tvar[4]!=1)*(j+2)+(tvar[4]==1)]) %*% solve(t(Tt[,,(tvar[2]!=1)*(j+2)+(tvar[2]==1)]))
        if(j==ss$n){break}else{
            # Calculating log likelihood value
            loglik <- loglik + dmvnorm(x = ss$y[j+1,], mean = matrix( Zt[,,(tvar[1]!=1)*(j+1)+(tvar[1]==1)], nrow = 1) %*% Tt[,,(tvar[2]!=1)*(j+1)+(tvar[2]==1)] %*% matrix( ss$mf[j,], ncol = 1), sigma = matrix( Zt[,,(tvar[1]!=1)*(j+1)+(tvar[1]==1)], nrow = 1) %*% (Tt[,,(tvar[2]!=1)*(j+1)+(tvar[2]==1)] %*% ss$Cf[[j]] %*% t(Tt[,,(tvar[2]!=1)*(j+1)+(tvar[2]==1)]) + Qt[,,(tvar[4]!=1)*(j+1)+(tvar[4]==1)]) %*% matrix( Zt[,,(tvar[1]!=1)*(j+1)+(tvar[1]==1)], ncol = 1) + Ht[,,(tvar[3]!=1)*(j+1)+(tvar[3]==1)], log = TRUE  )
        }
    }

        ss$m <- t(kfas.filt$alphahat)
        ss$C <- list()
        for(j in 1:ss$n){
            ss$C[[j]] <- kfas.filt$V[,,j]
        }

    ss$loglik <- loglik

    list(kfas=kfas.filt,ss=ss)
}

sspir2kfas <- function(M, dim, ...){
    ## Tiny help function
    Mt <- vapply(1:dim[3], M, array(0,dim=c(dim[1],dim[2])), ...)
}



Fkfs <- function(ss, tvar, offset=1) {
  if(ss$family$family=="gaussian"){
    return(Fkfilter(ss, tvar))
  }else{
      ## Check of time-varying matrices:
      Zt <- sspir2kfas(ss$Fmat, dim=c(ss$d, ss$p, tvar[1]), x=ss$x, phi=ss$phi)
      Tt <- sspir2kfas(ss$Gmat, dim=c(ss$p, ss$p, tvar[2]), x=ss$x, phi=ss$phi)
      Qt <- sspir2kfas(ss$Wmat, dim=c(ss$p, ss$p, tvar[4]), x=ss$x, phi=ss$phi)

      # Changing initial values to fit to KFAS
      a1 <- Tt[,,1] %*% matrix(ss$m0, ncol=1)
      P1 <- (Tt[,,1] %*% ss$C0 %*% t(Tt[,,1]) + Qt[,,1])

      # Changing letter in family to fit to KFAS notation
      if(ss$family$family == "poisson"){
          dist <- "Poisson"
      }else{
          if(ss$family$family == "binomial"){
              dist <- "Binomial"
          }else{
              stop("Family not allowed!")
          }
      }

      # specify model in KFAS notation
      m1 <- KFAS::SSModel(y=ss$y, Z=Zt, T=Tt, R=diag(ss$p), Q=Qt, a1=a1, P1=P1, u=offset, distribution=dist)
      # Calculates approximating Gaussian model
      out <- KFAS::approxSSM(m1)
      # specify approximating Gaussian model
      mg <- m1
      mg$y <- ss$ytilde <- out$y
      mg$H <- ss$vtilde <- out$H
      mg$distribution <- "Gaussian"
      outg <- KFAS::KFS(mg, smoothing="state", simplify=FALSE)

      # Backtracking filtered means, m, and covariances, C
      ss$mf <- t( (solve(Tt[,,1:tvar[2]]) %*% outg$a)[,-c(1)] )
      ss$Cf <- list()      
      for(j in 1:ss$n){
          ss$Cf[[j]] <- solve(Tt[,,(tvar[2]!=1)*(j+2)+(tvar[2]==1)]) %*% (outg$P[,,j+1] - Qt[,,(tvar[4]!=1)*(j+2)+(tvar[4]==1)]) %*% solve(t(Tt[,,(tvar[2]!=1)*(j+2)+(tvar[2]==1)]))
      }

      # Backtracking smoothed means, m.tilde, and covariances, C.tilde
      ss$m <- t(outg$alphahat)
      ss$C <- list()
      for(j in 1:ss$n){
          ss$C[[j]] <- outg$V[,,j]
      }

      ss$mu     <- out$theta

      ss$loglik <- outg$logLik

      out$app.model   <- outg

      list(kfas=out, ss=ss)
  }
}


######################################################################

"smoother" <-
function(ss) {
  UseMethod("smoother")
}

"smootherstep" <-
function(m,C,Gmatx,Wtx,mx,Cx)
  {
    Rx <- Gmatx %*% C %*% t(Gmatx) + Wtx# Gmatx = G_{t+1}
#    print(Rx)
Rx <- (Rx + t(Rx))/2
    B  <- C %*% t(Gmatx) %*% solve( Rx )
    ms <- t(m) + B%*%(t(mx) - Gmatx %*% t(m)) # mx=ms_{t+1}
    Cs <- C + B%*%(Cx-Rx)%*%t(B)
    assign(".Last.mtilde", m, envir=.GlobalEnv)
    list(ms=ms,Cs=Cs)
  }

"smootherstep.uni" <-
function(m,C,Gmatx,Wtx,mx,Cx)
  {
    Rx <- Gmatx * C * Gmatx + Wtx# Gmatx = G_{t+1}
    B  <- C * Gmatx / Rx
    ms <- m + B*(mx - Gmatx * m) # mx=ms_{t+1}
    Cs <- C + B*(Cx-Rx)*B
    assign(".Last.mtilde", m, envir=.GlobalEnv)
    list(ms=ms,Cs=Cs)
  }



"smoother.SS" <-
function(ss) {

  mf <- ss$m
  Cf <- ss$C
  m0f<- ss$m0
  C0f<- ss$C0

  m <- ss$m
  C <- ss$C
  m0<- ss$m0
  C0<- ss$C0

  d <- ss$d

  ## The state space model is given as
  ## Observation eq:    y[,t] = Fmat^T * theta_t  + N_d(0,diag(V[,t]))
  ## Evolution   eq:  tht_t   = Gmat * theta_{t-1}+ N_q(0,diag(W[,t]))
  ## Initial       :  tht_0   ~ N(m0,C0)

  ## The Kalman filter yields m_t and C_t, with
  ##            (tht_t|y[,1:t]) ~ N(m_t, C_t)

  ## The Kalman smoother yields ms_t and Cs_t, with
  ##            (tht_t|y[,1:nobs]) ~ N(ms_t, Cs_t)

  ##
  ## In the general case, the (backwards) recursion goes:
  ##  R_{t+1}= G_{t+1} * C_t * G_{t+1}^T + W_{t+1}
  ##  B_t = C_t * G_{t+1}^T * R_{t+1}^{-1}
  ##
  ##  ms_t = m_t + B_t * ( ms_{t+1} - G_t+1 * m_t)
  ##  Cs_t = C_t + B_t * ( Cs_{t+1} - (R_t+1) )*B_t^T

  ##
  ## in our case, F_t = Fmat, G_t = Gmat, V_t = diag(V[,t]),
  ##              W_t = diag(W[,t])

  ## y: matrix with T columns each holding a d-dimensional vector
  ## Gmat: qxq matrix
  ## Wt:   qxT matrix holding variances of tht. in columns
  ## m: qxnobs   from Kalman filter
  ## C: qxqxnobs from Kalman filter

  nobs <- ss$n
#  q    <- dim(m)[1]

    ## Check consistency in input
#    if (dim(C)[1] != q || dim(C)[2] != q || dim(C)[3] != nobs)
#        stop("Dimension of Fmat (",dim(Fmat)[1]," by ", dim(Fmat)[2],
#            ") should be ",q," by ",d,"!\n")
#    if (dim(Gmat)[1] != q || dim(Gmat)[2] != q)
#        stop("Dimension of Gmat (",dim(Gmat)[1]," by ", dim(Gmat)[2],
#            ") should be ",q," by ",q,"!\n")
#    if (dim(Wt)[1] != q || dim(Wt)[2] != nobs)
#       stop("Dimension of Wt (",dim(Wt)[1]," by ", dim(Wt)[2],
#          ") should be ",q," by ",nobs,"!\n")

    ## output matrices are just m and C, since the backwards recursion
    ## only uses 'future' m's and C's

        ## Nothing to do in step 1 (ms[,nobs] = m[,nobs])

    ## run the recursion

  mu <- matrix(NA, nobs, d)
   mu[nobs,] <- t(ss$Fmat(nobs,ss$x,ss$phi)) %*% m[nobs,]

  for (tt in (nobs-1):1)
    {
#      cat(tt," ",sep="")
if(class(ss$Wmat)=="matrix")
	{
      if (ss$p == 1)
        nextstep <-
        smootherstep.uni(
                     m[tt,],
                     C[[tt]],
                     ss$Gmat(tt+1,ss$x,ss$phi),
                     ss$Wmat,
                     m[tt+1,],
                     C[[tt+1]])

      else
      nextstep <-
        smootherstep(
                     matrix(m[tt,],nrow=1),
                     C[[tt]],
                     ss$Gmat(tt,ss$x,ss$phi),
                     ss$Wmat,
                     matrix(m[tt+1,],nrow=1),
                     C[[tt+1]])
	}
	else {
      if (ss$p == 1)
        nextstep <-
        smootherstep.uni(
                     m[tt,],
                     C[[tt]],
                     ss$Gmat(tt+1,ss$x,ss$phi),
                     ss$Wmat(tt+1,ss$x,ss$phi),
                     m[tt+1,],
                     C[[tt+1]])

      else
      nextstep <-
        smootherstep(
                     matrix(m[tt,],nrow=1),
                     C[[tt]],
                     ss$Gmat(tt,ss$x,ss$phi),
                     ss$Wmat(tt+1,ss$x,ss$phi),
                     matrix(m[tt+1,],nrow=1),
                     C[[tt+1]])
	}

      m[tt,]  <- nextstep$ms
      C[[tt]] <- nextstep$Cs
      mu[tt,] <- t(ss$Fmat(tt,ss$x,ss$phi)) %*% m[tt,]
      }
    ## do the last step
#  if (ss$p == 1)
#    laststep <- smootherstep.uni(
#                             m0,
#                             C0,
#                             ss$Gmat(1,ss$x,ss$phi),
#                             ss$Wmat(1,ss$x,ss$phi),
#                             m[,1],
#                             C[[1]]
#                             )
#  else
#    laststep <- smootherstep(
#                             m0,
#                             C0,
#                             ss$Gmat(1,ss$x,ss$phi),
#                             ss$Wmat(1,ss$x,ss$phi),
#                             m[,1],
#                             C[[1]]
#                             )
#    m0 <- laststep$ms
#    C0 <- laststep$Cs

  ss$m0 <- m0
  ss$C0 <- C0
  if (is.ts(ss$y))
    ss$m <- ts(m,start(ss$y),end=end(ss$y),frequency=frequency(ss$y))
  else
    ss$m <- m
  ss$C <- C
  ss$mu <- mu

  ss$mf <- mf
  ss$Cf <- Cf
  ss$m0f<- m0f
  ss$C0f<- C0f

  class(ss) <- c("Smoothed","SS")
  ss
  }

smoother.ssm <- function(ss) {
  return(kfs(ss))
}
