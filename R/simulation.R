"ksimulate" <-
function(ss,N=1) {
  UseMethod("ksimulate")
}

"ksimulate.SS" <-
function(ss,N=1) {
#  require(MASS)
  ## Assumes m,C to be from the FILTER (not the smoother!)
  
  ## calculates R_t=G_t*C_{t-1}*G_t^t + W_t,
  ## a_t = G_t m_{t-1} and B_t = C_t G_{t+1}^t R_{t+1}^-1
  ## from m_t, C_t, Gmat, Wmat. 

  ## Then simulates from theta_n |y ~ N(m_n, C_n) and downwards recursion
  ## theta_t | theta_{t+1}  ~ N(m_t+B_t(theta_{t+1}-a_{t+1}, C_t - B_tR_{t+1}B_t^t)

  ## N simulations are performed and are returned in a list
  ## simul[[1]]..simul[[N]] which are concatenated to the ss-object.

  n <- ss$n
  p <- ss$p
  Gmat <- ss$Gmat
  Wmat <- ss$Wmat
  m    <- ss$m
  m0   <- ss$m0
  C    <- ss$C
  C0    <- ss$C0
  
  R <- vector("list",n)
  B <- vector("list",n)
#  a <- matrix(NA,p,n)
  a <- matrix(NA,n,p)
  
  GG <- Gmat(1, ss$x, ss$phi)

  if (ss$p == 1)
    {
      cat("Univariate\n")
      a[1,] <- GG * m0
      R[[1]]<- GG * C0 * GG + Wmat(1,ss$x,ss$phi)
    }
  else
    {
      a[1,] <- GG %*% t(m0)
      R[[1]]<- GG %*% C0 %*% t(GG) + Wmat(1,ss$x,ss$phi)
    }
  for (tt in 2:n) {
    ## calculate R, a
    GG <- Gmat(tt,ss$x,ss$phi)
    if (ss$p == 1)
      {
        a[tt,] <- GG * m[tt-1,]
        R[[tt]]<- GG * C[[tt-1]] * GG + Wmat(tt,ss$x,ss$phi)
        class(R[[tt]]) <- c("Hermitian","Matrix") # p.d.
  }
    else
      {
        a[tt,] <- GG %*% matrix(m[tt-1,])
        R[[tt]]<- GG %*% C[[tt-1]] %*% t(GG) + Wmat(tt,ss$x,ss$phi)
        class(R[[tt]]) <- c("Hermitian","Matrix") # p.d.
      }
  }

  simul <- array(NA,dim=c(n,p,N))
  ## simul at time n

  if (ss$p == 1) 
    simul[n,,] <- rnorm(N, m[n,], sqrt(C[[n]]))
  else 
    simul[n,,] <- t(mvrnorm(N, m[n,], C[[n]]))
  
  for (tt in (n-1):1) {
    ## simulate at time tt

    if (ss$p == 1)
      {
        Bt <- C[[tt]] * Gmat(tt+1,ss$x,ss$phi) / R[[tt+1]]
        Vart<-C[[tt]] - Bt*R[[tt+1]]*Bt
        for (it in 1:N) 
          simul[tt,,it] <- rnorm(1,
                                   m[tt,] + Bt*(simul[tt+1,,it]-a[tt+1,]),
                                   sqrt(Vart)
                                   )
      }
    else
      {
        Bt <- C[[tt]] %*% t(Gmat(tt+1,ss$x,ss$phi)) %*% solve(R[[tt+1]])
        Vart<-C[[tt]] - Bt%*%R[[tt+1]]%*%t(Bt)
        for (it in 1:N) 
          simul[tt,,it] <- t(mvrnorm(1,
                                   m[tt,] + Bt%*%(simul[tt+1,,it]-a[tt+1,]),
                                   Vart
                                   ))
      }
  }
  
  simul
}
