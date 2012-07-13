\name{Fkfs}
\alias{Fkfs}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fast Kalman filtering and smoothing of state space models
  using \code{\link[KFAS]{KFS}} or \code{\link[KFAS]{logLik.SSModel}}}
\description{
  Does Kalman filtering followed by smoothing using
  \code{\link[KFAS]{KFS}} for gaussian state
  space models otherwise by use of \code{\link[KFAS]{logLik.SSModel}}.
}
\usage{
Fkfs(ss, tvar, offset = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{object of class \code{\link{ssm}}.}
  \item{tvar}{vector with four elements indicating the dynamic nature of
  matrices in the object \code{ssm}. The first element corresponds to
  \code{ss$Fmat}, the second to \code{ss$Gmat}, the third to
  \code{ss$Vmat} and the fourth to \code{ss$Wmat}. If \code{tvar[i]==1}
  the corresponding matrix is constant, otherwise \code{tvar[i]} is set
  to \code{ss$n}.}
  \item{offset}{vector of length equal to length of observations.}
}
\details{
The Kalman filter and smoother in KFAS provides faster computations than \code{\link{kfilter}} and \code{\link{smoother}}, however requires alternative specification of the state space model. This function transform the object \code{SS}on which Kalman filtering followed by smoothing using \code{\link[KFAS]{KFS}} may be performed.
}
\value{List with two objects. The first object, \code{kfas}, is the output
  from \code{\link[KFAS]{KFS}}. The second object, \code{ss}, of class \code{\link{SS}} with updated \code{m}, \code{C} and \code{loglik}.
}
\author{Anette Luther Christensen and Claus Dethlefsen.}
\seealso{\code{\link{Fkfilter}}}
\examples{
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
}
\keyword{models}% at least one, from doc/KEYWORDS
