pkgname <- "sspir"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('sspir')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("EMalgo")
### * EMalgo

flush(stderr()); flush(stdout())

### Name: EMalgo
### Title: Estimation of variance matrices in a Gaussian state space model
### Aliases: EMalgo
### Keywords: models

### ** Examples

## Formulates Gaussian state space model:
## Trend: local linear
## Seasonal variation: sine function with frequency = 1
m1 <- SS( Fmat = function(tt, x, phi) {
            Fmat      <- matrix(NA, nrow=3, ncol=1)
            Fmat[1,1] <- 1
            Fmat[2,1] <- cos(2*pi*tt/12)
            Fmat[3,1] <- sin(2*pi*tt/12)
            return(Fmat)
          },
          Gmat = function(tt, x, phi) {
            matrix(c(1,0,0,0,1,0,0,0,1), nrow=3)
          },
          Wmat = matrix(c(0.01,0,0,0,0.1,0,0,0,0.1), nrow=3),
          Vmat = matrix(1),
          m0   = matrix(c(0,0,0), nrow=1),
          C0   = matrix(c(1,0,0,0,0.001,0,0,0,0.001), nrow=3, ncol=3)
        )

## Simulates 100 observation from m1
m1 <- recursion(m1, 100)

## Specifies the correlation structure of W:
Wstruc <- function(W){	
            W[1,2:3] <- 0
            W[2:3,1] <- 0
            W[2,2]   <- (W[2,2]+W[3,3])/2
            W[3,3]   <- W[2,2]
            W[2,3]   <- 0
            W[3,2]   <- 0
            return(W)
          }

## Estimstes variances and covariances of W by use of the EM algorithm
estimates <- EMalgo(m1, Wstruc=Wstruc)
estimates$Wmat.est

## Plots estimated model
plot(estimates$model)



cleanEx()
nameEx("Fkfilter")
### * Fkfilter

flush(stderr()); flush(stdout())

### Name: Fkfilter
### Title: Kalman filtering of Gaussian state space model using'KFS'
### Aliases: Fkfilter
### Keywords: models

### ** Examples

## Formulates Gaussian state space model:
## Trend: local linear
## Seasonal variation: sine function with frequency = 1
m1 <- SS( Fmat = function(tt, x, phi) {
            Fmat      <- matrix(NA, nrow=3, ncol=1)
            Fmat[1,1] <- 1
            Fmat[2,1] <- cos(2*pi*tt/12)
            Fmat[3,1] <- sin(2*pi*tt/12)
            return(Fmat)
          },
          Gmat = function(tt, x, phi) {
            matrix(c(1,0,0,0,1,0,0,0,1), nrow=3)
          },
          Wmat = function(tt, x, phi) {matrix(c(0.01,0,0,0,0.1,0,0,0,0.1), nrow=3)},
          Vmat = function(tt, x, phi) {matrix(1)},
          m0   = matrix(c(0,0,0), nrow=1),
          C0   = matrix(c(1,0,0,0,0.001,0,0,0,0.001), nrow=3, ncol=3)
        )

## Simulates 100 observation from m1
m1 <- recursion(m1, 100)

filtered <- Fkfilter(m1, tvar=c(m1$n, 1, 1, 1))

plot(filtered$ss)




cleanEx()
nameEx("Fkfs")
### * Fkfs

flush(stderr()); flush(stdout())

### Name: Fkfs
### Title: Fast Kalman filtering and smoothing of state space models using
###   'KFS' or 'logLik.SSModel'
### Aliases: Fkfs
### Keywords: models

### ** Examples

## Gaussian model to simulate observations from
m1 <- SS( Fmat = function(tt, x, phi) {
            Fmat      <- matrix(NA, nrow=3, ncol=1)
            Fmat[1,1] <- 1
            Fmat[2,1] <- cos(2*pi*tt/12)
            Fmat[3,1] <- sin(2*pi*tt/12)
            return(Fmat)
          },
          Gmat = function(tt, x, phi) {
            matrix(c(1,0,0,0,1,0,0,0,1), nrow=3)
          },
          Wmat = function(tt, x, phi) {
                    matrix(c(0.01,0,0,0,0.1,0,0,0,0.1), nrow=3)},
          Vmat = function(tt, x, phi) {matrix(1)},
          m0   = matrix(c(0,0,0), nrow=1),
          C0   = matrix(c(1,0,0,0,0.001,0,0,0,0.001), nrow=3, ncol=3)
        )

## Simulates 100 observation from m1
m1 <- recursion(m1, 100)

## Formulates model og class ssm
ssm1 <- ssm(m1$y ~ -1 + tvar(polytime(1:m1$n),1)
                      + tvar(polytrig(1:m1$n, 1)), 
                        family="gaussian", fit=FALSE)

smoothed <- Fkfs(ssm1$ss, tvar=c(m1$n, 1, 1, 1))

## Non Gaussian model
phi.start <- StructTS(log10(UKgas),type="BSM")$coef[c(4,1,2,3)]
gasmodel <- ssm( log10(UKgas) ~ -1+ tvar(polytime(time,1))+
tvar(sumseason(time,4)), phi=phi.start, fit=FALSE)

smoothed <- Fkfs(gasmodel$ss,tvar=c(1,1,1,1))
## Trend plot
## Trend plot
ts.plot( smoothed$kfas$alphahat[1,] )
## Season plot
ts.plot( smoothed$kfas$alphahat[3,] )



cleanEx()
nameEx("SS")
### * SS

flush(stderr()); flush(stdout())

### Name: SS
### Title: Representation of Gaussian State Space Model
### Aliases: SS plot.SS print.SS SS-class C0.SS C0<-.SS Fmat.SS Fmat<-.SS
###   Gmat.SS Gmat<-.SS m0.SS m0<-.SS phi.SS phi<-.SS Vmat.SS Vmat<-.SS
###   Wmat.SS Wmat<-.SS
### Keywords: models

### ** Examples

data(kurit)  ## West & Harrison, page 40
m1 <- SS(y=kurit,
         Fmat=function(tt,x,phi) return(matrix(1)),
         Gmat=function(tt,x,phi) return(matrix(1)),
         Wmat=function(tt,x,phi) return(matrix(5)), ## Alternatively Wmat=matrix(5)
         Vmat=function(tt,x,phi) return(matrix(100)), ## Alternatively Vmat=matrix(100)
         m0=matrix(130),C0=matrix(400)
         )

plot(m1$y)
m1.f <- kfilter(m1)
m1.s <- smoother(m1.f)
lines(m1.f$m,lty=2,col=2)
lines(m1.s$m,lty=2,col=2)

## make a model with an intervention at time 10
m2 <- m1
Wmat(m2) <- function(tt,x,phi) {
  if (tt==10) return(matrix(900))
  else return(matrix(5))
}

m2.f <- kfilter(m2)
m2.s <- smoother(m2.f)
lines(m2.f$m,lty=2,col=4)
lines(m2.s$m,lty=2,col=4)

## Use 'ssm' to construct an SS skeleton
phi.start <- StructTS(log10(UKgas),type="BSM")$coef[c(4,1,2,3)]
gasmodel <- ssm( log10(UKgas) ~ -1+
                 tvar(polytime(time,1))+
                 tvar(sumseason(time,12)),
                 phi=phi.start)

m0(gasmodel)
C0(gasmodel)
phi(gasmodel)

fit <- getFit(gasmodel)
plot( fit$m[,1:3]  )



cleanEx()
nameEx("extended")
### * extended

flush(stderr()); flush(stdout())

### Name: extended
### Title: Iterated Extended Kalman Smoothing
### Aliases: extended extended.SS print.extended
### Keywords: models

### ** Examples

data(mumps)
index <- 1:length(mumps) # use 'index' instead of time
model <- ssm( mumps ~ -1 + tvar(polytime(index,1)),
              family=poisson(link=log))
results <- getFit(model)
plot(mumps,type='l',ylab='Number of Cases',xlab='',axes=FALSE)
lines( exp(results$m[,1]), lwd=2)
## Alternatives:
## results2 <- extended(model$ss)
## results3 <- kfs(model) ## yields the same



cleanEx()
nameEx("kfilter")
### * kfilter

flush(stderr()); flush(stdout())

### Name: kfilter
### Title: Kalman filter for Gaussian state space model
### Aliases: kfilter kfilter.SS filterstep
### Keywords: models

### ** Examples

data(kurit)
m1 <- SS(kurit)
phi(m1) <- c(100,5)
m0(m1) <- matrix(130)
C0(m1) <- matrix(400)

m1.f <- kfilter(m1)
plot(m1$y)
lines(m1.f$m,lty=2)



cleanEx()
nameEx("kfs")
### * kfs

flush(stderr()); flush(stdout())

### Name: kfs
### Title: (Iterated extended) Kalman smoother
### Aliases: kfs kfs.ssm kfs.SS
### Keywords: models

### ** Examples

data(mumps)
index <- 1:length(mumps)
phi.start <- c(0,0,0.0005,0.0001)
m3 <- ssm( mumps ~ -1 + tvar(polytime(index,1)) +
                  tvar(polytrig(index,12,1)),
                  family=poisson(link=log),
                  phi=phi.start, C0 = diag(4),
                  fit=FALSE
)

## The option "fit=FALSE" means that the Kalman Filter/Smoother is not
## run.
## At this point you may inspect/change the setup before running 'kfs'
C0(m3)
C0(m3) <- 10*diag(4)
## incorporate possible structural 'jump' at timepoint 10
Wold <- Wmat(m3)
Wmat(m3) <- function(tt,x,phi) {
    W <- Wold(tt,x,phi)
    if (tt==10) {W[2,2] <- 100*W[2,2]; return(W)}
    else return(W)
}

m3.fit <- kfs(m3)

plot(mumps,type='l',ylab='Number of Cases',xlab='')
lines(exp(m3.fit$m[,1]),type='l',lwd=2)



cleanEx()
nameEx("recursion")
### * recursion

flush(stderr()); flush(stdout())

### Name: recursion
### Title: Simulate from a Gaussian state space model
### Aliases: recursion recursion.SS
### Keywords: models

### ** Examples

data(kurit)
m1 <- SS(kurit)
phi(m1) <- c(100,5)
m0(m1) <- matrix(130)
C0(m1) <- matrix(400)

par(mfrow=c(2,1))
plot(recursion(m1,100))
phi(m1) <- c(5,100)
plot(recursion(m1,100))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("simulate")
### * simulate

flush(stderr()); flush(stdout())

### Name: ksimulate
### Title: Forwards filtering Backwards sampling
### Aliases: ksimulate ksimulate.SS
### Keywords: models

### ** Examples

data(kurit)
m1 <- SS(kurit)
phi(m1) <- c(100,5)
m0(m1) <- matrix(130)
C0(m1) <- matrix(400)

m1 <- kfilter(m1)
m1.s <- smoother(m1)
sim <- ksimulate(m1,10)

plot(kurit)
for (i in 1:10) lines(sim[,,i],lty=2,col=2)

lines(smoother(m1)$m,lwd=2)



cleanEx()
nameEx("smoother")
### * smoother

flush(stderr()); flush(stdout())

### Name: smoother
### Title: Kalman smoother for Gaussian state space model
### Aliases: smoother smoother.SS smootherstep smootherstep.uni
###   print.Smoothed
### Keywords: models

### ** Examples

data(kurit)
m1 <- SS(kurit)
phi(m1) <- c(100,5)
m0(m1) <- matrix(130)
C0(m1) <- matrix(400)

m1.s <- smoother(kfilter(m1))
plot(m1$y)
lines(m1.s$m,lty=2)



cleanEx()
nameEx("specials")
### * specials

flush(stderr()); flush(stdout())

### Name: Specials
### Title: Special functions used in ssm formulas
### Aliases: polytrig polytime sumseason season
### Keywords: models

### ** Examples

polytrig(1:10,degree=2)
season(1:12,period=4)



cleanEx()
nameEx("ssm")
### * ssm

flush(stderr()); flush(stdout())

### Name: ssm
### Title: Define state-space model in a glm-style call.
### Aliases: ssm kfilter.ssm smoother.ssm C0 C0<- Fmat Fmat<- Gmat Gmat<-
###   m0 m0<- phi phi<- Vmat Vmat<- Wmat Wmat<- C0.ssm C0<-.ssm Fmat.ssm
###   Fmat<-.ssm Gmat.ssm Gmat<-.ssm m0.ssm m0<-.ssm phi.ssm phi<-.ssm
###   Vmat.ssm Vmat<-.ssm Wmat.ssm Wmat<-.ssm getFit print.ssm
### Keywords: models

### ** Examples

data(vandrivers)
vandrivers$y <- ts(vandrivers$y,start=1969,frequency=12)
vd.time <- time(vandrivers$y)
vd <- ssm( y ~ tvar(1) + seatbelt + sumseason(vd.time,12),
          family=poisson(link="log"),
          data=vandrivers,
          phi = c(1,0.0004),
          C0=diag(13)*100,
          fit=FALSE
          )
phi(vd)["(Intercept)"] <- exp(- 2*3.703307 )
C0(vd) <- diag(13)*1000
vd.res <- kfs(vd)

plot( vd.res$m[,1:3] )

attach(vandrivers)
plot(y,ylim=c(0,20))
lines(exp(vd.res$m[,1]+vd.res$m[,2]*seatbelt),lwd=2 )
detach(vandrivers)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
