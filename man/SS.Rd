\name{SS}
\alias{SS}
\alias{plot.SS}
\alias{print.SS}
\alias{SS-class}
\alias{C0.SS}
\alias{C0<-.SS}
\alias{Fmat.SS}
\alias{Fmat<-.SS}
\alias{Gmat.SS}
\alias{Gmat<-.SS}
\alias{m0.SS}
\alias{m0<-.SS}
\alias{phi.SS}
\alias{phi<-.SS}
\alias{Vmat.SS}
\alias{Vmat<-.SS}
\alias{Wmat.SS}
\alias{Wmat<-.SS}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Representation of Gaussian State Space Model}
\description{
  Creates an SS-object describing a Gaussian state space model.
}
\usage{
SS(y = NA, x = list(x = NA),
   Fmat = function(tt,x,phi) { return(matrix(1)) },
   Gmat = function(tt,x,phi) { return(matrix(1)) },
   Vmat = function(tt,x,phi) { return(matrix(phi[1])) },
   Wmat = function(tt,x,phi) { return(matrix(phi[2])) },
   m0 = matrix(0),
   C0 = matrix(100),
   phi = c(1,1))

\method{C0}{SS}(ssm) 
\method{m0}{SS}(ssm) 
\method{Fmat}{SS}(ssm)
\method{Gmat}{SS}(ssm)
\method{Vmat}{SS}(ssm)
\method{Wmat}{SS}(ssm)
\method{phi}{SS}(ssm) 
\method{C0}{SS}(ssm) <- value
\method{m0}{SS}(ssm) <- value
\method{Fmat}{SS}(ssm) <- value
\method{Gmat}{SS}(ssm) <- value
\method{Vmat}{SS}(ssm) <- value
\method{Wmat}{SS}(ssm) <- value
\method{phi}{SS}(ssm) <- value

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a matrix giving a multivariate time series of
    observations. The observation at time \code{tt} is
    \code{y[tt,]}. The dimension of \code{y} is \eqn{n \times
      d}. Preferably, \code{y} is a \code{\link{ts}} object.}
  \item{x}{a list of entities (eg. covariates) passed as argument to the functions
    \code{Fmat}, \code{Gmat}, \code{Vmat}, and \code{Wmat}.}
  \item{Fmat}{a function depending on the parameter-vector \code{phi},
    covariates \code{x} and returns the \eqn{p \times d} design matrix at time
    \code{tt}.} 
  \item{Gmat}{a function depending on the parameter-vector \code{phi},
    covariates \code{x} and returns the \eqn{p \times p} evolution matrix at
    time \code{tt}.}
  \item{Vmat}{a function depending on the parameter-vector \code{phi},
    covariates \code{x} and returns the \eqn{d \times d} (positive definit)
    variance matrix at time \code{tt}. The matrix \code{Vmat} may be specified 
    as a constant matrix not depending on \code{phi} coded with \code{matrix}.}
  \item{Wmat}{a function depending on the parameter-vector \code{phi},
    covariates \code{x} and returns the \eqn{p \times p} (positive
    semidefinite) evolution variance matrix at time \code{tt}. The matrix \code{Wmat} may be specified 
    as a constant matrix not depending on \code{phi} coded with \code{matrix}.}
  \item{m0}{a \eqn{1 \times p} matrix giving the initial state.}
  \item{C0}{a \eqn{p \times p} variance matrix giving the variance matrix of
    the initial state.}
  \item{phi}{a (hyper) parameter vector passed as argument to the functions \code{Fmat}, \code{Gmat}, \code{Vmat}, and \code{Wmat}.}
  \item{ssm}{an object of class \code{SS}.}
  	\item{value}{an object to be assigned to the element of the
	  state space model.}

}
\details{
    The state space model is given by
    \deqn{Y_t     = F_t^\top  \theta_t     + v_t, \quad v_t \sim N(0,V_t)}{Y_t     = F_t^T * \theta_t     + v_t, v_t ~ N(0,V_t)}
    \deqn{\theta_t = G_t   \theta_{t-1} + w_t, \quad w_t \sim N(0,W_t)}{\theta_t = G_t  * \theta_{t-1} + w_t, w_t ~ N(0,W_t)}
    for \eqn{t=1,\ldots,n}{t=1,...,n}. The matrices \eqn{F_t}, \eqn{G_t}, \eqn{V_t},
    and \eqn{W_t} may depend on a parameter vector \eqn{\phi}.
    The initialization is given as
\deqn{\theta_0 \sim N(m_0,C_0).}{\theta_0 ~ N(m_0,C_0).}
}
\value{
  An object of class \code{SS}, which is a list with the following components
  \item{y}{as input.}
  \item{x}{as input.}
  \item{Fmat}{as input.}
  \item{Gmat}{as input.}
  \item{Vmat}{as input.}
  \item{Wmat}{as input.}
  \item{m0}{as input.}
  \item{C0}{as input.}
  \item{phi}{as input.}
  \item{n}{the number of time points}
  \item{d}{the dimension of each observation.}
  \item{p}{the dimension of the state vector at each timepoint.}
  \item{ytilde}{adjusted observations for use in the extended Kalman
    filter, see \code{\link{extended}}.}
  \item{iteration}{an integer giving the number of iterations used in
    the extended Kalman filter, see \code{\link{extended}}.} 
  \item{m}{after Kalman filtering (or smoothing), holds the conditional
    mean of the state vectors given the observations up till time \eqn{t}
    (filtering) or all observations (smoothing). This is organised in a
    \eqn{n \times p} dimensional matrix holding \eqn{m_t} (\eqn{m_t^*})
    in rows. Is returned as a \code{\link{ts}} object.}
  \item{C}{after Kalman filtering (or smoothing), holds the conditional
    variance of the state vectors given the observations up til time \eqn{t}
    (filtering) or all observations (smoothing). This is organised in a
    list holding the \eqn{p \times p} dimensional matrices \eqn{C_t}
  (\eqn{C_t^*}).} 
  \item{mu}{after Kalman smoothing, holds the conditional
    mean of the signal (\eqn{\mu_t=F_t^\top \theta_t}) given all
    observations. This is organised in a \eqn{n \times d} dimensional
    matrix holding \eqn{\mu_t} in rows.} 
  \item{loglik}{the log-likelihood value after Kalman filtering.}
}

\author{Claus Dethlefsen, \enc{Søren}{Soren} Lundbye-Christensen and Anette Luther Christensen}

\seealso{\code{\link{ssm}} for a glm-like interface of specifying
  models, \code{\link{kfilter}} for Kalman filter and
  \code{\link{smoother}} for Kalman smoother.}
\examples{
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
}
\keyword{models}% at least one, from doc/KEYWORDS
