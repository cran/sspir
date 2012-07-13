\name{sspir-package}
\alias{sspir-package}
\alias{sspir}
\docType{package}
\title{State Space Models in R}
\description{The main contribution of this package is to give a formula
  language for specifying dynamic generalized linear models. That is, an
  extension of glm formulae by marking terms with \code{tvar} to specify
  that their coefficients are time-varying. The package also provides
  (extended) Kalman filter and Kalman smoother for models within
  Gaussian, Poisson and binomial families. To get started, try
  \code{demo(gas)}, \code{demo(vandrivers)} and \code{demo(mumps)}.

  Also, read the paper "C. Dethlefsen and
  S. Lundbye-Christensen. Formulating state space models in r with focus
  on longitudinal regression models. Journal of Statistical Software,
  16(1), 2006." 
}

\section{What sspir does not include}{
  \itemize{
  \item{Optimization}{
    To keep full generality, the Kalman filter and
    smoother use the \code{\link{SS}} functions \code{Fmat},
    \code{Gmat}, \code{Vmat}, \code{Wmat}. The special cases where these
    matrices are time-invariant, the algorithms can be considerably
    speeded up by implementing these special cases eg. in C. For simple
    models, see \code{\link[stats]{StructTS}}.}
    \item{Diffuse initialization}{ 
    We use \eqn{m_0} and \eqn{C_0} as
    initialization of the state process. This may cause numeric problems
  which may be solved using diffuse initialization (see Durbin and
  Koopman (2000)).}
\item{ARIMA models}{
  These are not directly supported, but since they
  can be expressed as state space models, it is possible to specify
  them as \code{\link{SS}} objects.}
\item{Importance sampling}{
  We use iterated extended Kalman smoothing
  and use the likelihood from the approximating Gaussian state space
  model. This may be improved using importance sampling.}
\item{Multivariate observations}{
  We have plans including methods for
  combining \code{\link{ssm}} objects for different time-series so
  that they can possibly share components in the latent process.}
\item{Multi-process models}{
  We have plans for including
  methods for combining \code{\link{ssm}} objects defining different
  models for the same time-series. Providing prior probabilities for
  these models, it is possible to calculate posterior probabilites for
  the models, thereby discriminating between the models.}
\item{Markov chain Monte Carlo}{
  A state space model may be part of a
  hierarchical model giving priors on hyper-parameters. Inference may
  be done using MCMC methods and for the state space model, the
  Forwards filtering, Backwards sampling method may be used. We plan
  to include at least some support for this.}
\item{Missing values}{
  Missing values are not allowed in the covariates.}
\item{Family}{Currently, only three combinations of distribution/link
  functions are supported: Gaussian/identity, Poisson/log,
  Binomial/logit. To add new variations, edit the function \code{getFamily}}
}
}

\section{Other software}{
  \itemize{
  \item{Ox/SsfPack}{
    SsfPack is a package for Ox by SJ Koopman,
    N Shephard and JA Doornik. Multivariate Gaussian state
    space models with some support for non-Gaussian
    models. \url{http://www.ssfpack.com/}. Closely related to
    S+FinMetrics, a package for S-Plus, developed by Insightful. 
   \url{http://www.insightful.com/products/finmetrics/default.asp}.}
 \item{StructTS}{
   StructTS by BD Ripley is a part of R. Univariate
   Gaussian state space models of the class structural time series models. }
 \item{KEE}{
   Kalman Estimating Equations is a package for
   S-Plus by SJ Knudsen. Poisson-Tweedie log-linear
   State-Space
   model. \url{http://www.statdem.sdu.dk/~sjk/kee/index.html}.}
 \item{IEKS}{
   IEKS library is a package for S-Plus and \R by B Klein.
   Implementation of the iterated extended Kalman filter and smoother
   to discrete data from an exponential family. The latent process is
   assumed to be linear and
   Gaussian. \url{http://genetics.agrsci.dk/~bmk/}. }
 \item{Matlab}{
   Kalman filter toolbox for Matlab by K
   Murphy.
   \url{http://www.ai.mit.edu/~murphyk/Software/Kalman/kalman.html}. }
 \item{Bats}{
   BATS software by A Pole, M West and J
  Harrison.
  \url{http://www.isds.duke.edu/~mw/bats.html}. }
\item{dse}{
  Dynamic Systems Estimation library (dse1 and
   dse2) are packages for 
   R by P Gilbert. Multivariate Gaussian state space
   models with focus on ARMA
   models. See \url{http://www.bank-banque-canada.ca/pgilbert}.  }
 \item{CTSM}{
   Continuous Time Stochastic Modelling (CTSM) is a stand alone
   program by NR Kristensen, LE Christiansen and
   H Madsen. Stochastic differential equations, linear and
   non-linear Gaussian state space
   models. \url{http://www.imm.dtu.dk/ctsm/}.}
}
 Probably this list is incomplete. Feel free to contribute with
 more links to packages.
  }

\author{Claus Dethlefsen, \enc{Søren}{Soren} Lundbye-Christensen and Anette Luther Christensen.}
\seealso{\code{\link{ssm}}}
\keyword{models}% at least one, from doc/KEYWORDS
