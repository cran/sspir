\name{SS}
\alias{SS}
\alias{plot.SS}
\alias{print.SS}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Representation of Gaussian State Space Model}
\description{
  Creates an SS-object describing a Gaussian state space model.
}
\usage{
SS(y = NA, x = NA,
   Fmat = function(tt, x, phi) { NA },
   Gmat = function(tt, x, phi) { NA },
   Vmat = function(tt, x, phi) { NA },
   Wmat = function(tt, x, phi) { NA },
   m0 = 0, C0 = NA,
   phi = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{a matrix giving a multivariate time series of
    observations. The observation at time \code{tt} is
    \code{y[,tt]}. The dimension of \code{y} is \eqn{d \times n}.}
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
    variance matrix at time \code{tt}. }
  \item{Wmat}{a function depending on the parameter-vector \code{phi},
    covariates \code{x} and returns the \eqn{p \times p} (positive
    semidefinite) evolution variance matrix at time \code{tt}.}
  \item{m0}{a \eqn{p \times 1} matrix giving the initial state.}
  \item{C0}{a \eqn{p \times p} variance matrix giving the variance matrix of
    the initial state.}
  \item{phi}{a parameter vector passed as argument to the functions \code{Fmat}, \code{Gmat}, \code{Vmat}, and \code{Wmat}.}
}
\details{
    The state space model is given by
    \deqn{Y_t     = F_t^T  \theta_t     + v_t, v_t \sim N(0,V_t)}{Y_t     = F_t^T * \theta_t     + v_t, v_t ~ N(0,V_t)}
    \deqn{\theta_t = G_t   \theta_{t-1} + w_t, w_t \sim N(0,W_t)}{\theta_t = G_t  * \theta_{t-1} + w_t, w_t ~ N(0,W_t)}
    for \eqn{t=1,...,n}. The matrices \eqn{F_t}, \eqn{G_t}, \eqn{V_t},
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
  \item{ytilde}{for use in the extended Kalman filter.}
  \item{iteration}{for use in the extended Kalman filter.}
  \item{m}{after Kalman filtering (or smoothing), holds the conditional
    mean of the state vectors given the observations up til time \eqn{t}
    (filtering) or all observations (smoothing). This is organised in a
    \eqn{p \times n} dimensional matrix holding \eqn{m_t} (\eqn{m_t^*}) in columns.}
  \item{C}{after Kalman filtering (or smoothing), holds the conditional
    variance of the state vectors given the observations up til time \eqn{t}
    (filtering) or all observations (smoothing). This is organised in a
    list holding the \eqn{p \times p} dimensional matrices \eqn{C_t}
  (\eqn{C_t^*}).} 
  \item{mu}{after Kalman smoothing, holds the conditional
    mean of the signal (\eqn{\mu_t=F_t^\top \theta_t}) given all
    observations. This is organised in a \eqn{d \times n} dimensional
    matrix holding \eqn{\mu_t} in columns.} 
  \item{likelihood}{the log-likelihood value after Kalman filtering.}
}

\author{Claus Dethlefsen and S�ren Lundbye-Christensen}

\seealso{\code{\link{ssm}} for a glm-like interface of specifying
  models, \code{\link{kfilter}} for Kalman filter and
  \code{\link{smoother}} for Kalman smoother.}
\examples{
time <- 1:length(UKgas)
gasmodel <- ssm( log10(UKgas) ~ -1+
                 tvar(polytime(time,1))+
                 tvar(sumseason(time,12)),time=time)

gasmodel$ss$phi <- StructTS(log10(UKgas),type="BSM")$coef[c(4,1,2,3)]

fit <- kfs(gasmodel)

plot( ts( t(fit$m[1:3,]) ) )
}
\keyword{models}% at least one, from doc/KEYWORDS