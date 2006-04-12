\name{ssm}
\alias{ssm}
\alias{kfilter.ssm}
\alias{smoother.ssm}
\alias{C0}
\alias{C0<-}
\alias{Fmat}
\alias{Fmat<-}
\alias{Gmat}
\alias{Gmat<-}
\alias{m0}
\alias{m0<-}
\alias{phi}
\alias{phi<-}
\alias{Vmat}
\alias{Vmat<-}
\alias{Wmat}
\alias{Wmat<-}
\alias{C0.ssm}
\alias{C0<-.ssm}
\alias{Fmat.ssm}
\alias{Fmat<-.ssm}
\alias{Gmat.ssm}
\alias{Gmat<-.ssm}
\alias{m0.ssm}
\alias{m0<-.ssm}
\alias{phi.ssm}
\alias{phi<-.ssm}
\alias{Vmat.ssm}
\alias{Vmat<-.ssm}
\alias{Wmat.ssm}
\alias{Wmat<-.ssm}
\alias{getFit}
\alias{print.ssm}

\title{Define state-space model in a glm-style call.}
\description{
Use a glm-style formula and family arguments to setup a state space
model.
}
\usage{
ssm(formula, family = gaussian, data = list(), subset =
    NULL, fit = TRUE, phi = NULL, m0 = NULL, C0 = NULL,
                 Fmat = NULL, Gmat = NULL, Vmat = NULL, Wmat = NULL)
\method{C0}{ssm}(ssm) 
\method{m0}{ssm}(ssm) 
\method{Fmat}{ssm}(ssm)
\method{Gmat}{ssm}(ssm)
\method{Vmat}{ssm}(ssm)
\method{Wmat}{ssm}(ssm)
\method{phi}{ssm}(ssm) 
\method{C0}{ssm}(ssm) <- value
\method{m0}{ssm}(ssm) <- value
\method{Fmat}{ssm}(ssm) <- value
\method{Gmat}{ssm}(ssm) <- value
\method{Vmat}{ssm}(ssm) <- value
\method{Wmat}{ssm}(ssm) <- value
\method{phi}{ssm}(ssm) <- value
getFit(ssm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a formula with univariate response on the lefthand
    side. The righthand side is a sum of terms and the special functions
  \code{\link{sumseason}}, \code{\link{polytime}},
  \code{\link{polytrig}}, and \code{\link{season}} can be used. Terms
  can be marked by the \code{tvar}-function to create a term with
  time-varying coefficients. A special case is \code{tvar(1)} meaning a
  random walk.}
  \item{family}{a description of the error distribution and link function to
          be used in the model. This can be a character string naming a
          family function, a family function or the result of a call to
          a family function.  (See \code{\link{family}} for details of
	  family functions.)}  
  \item{data}{an optional data frame containing the variables in the model.
           If not found in \code{data}, the variables are taken from
          \code{environment(formula)}, typically the environment from which
          \code{ssm} is called.}
  \item{subset}{an optional vector specifying a subset of observations to be
          used in the fitting process.}
	\item{ssm}{an object of class \code{ssm}.}
	\item{fit}{a logical. If \code{TRUE}, the model is fit using the
	  \code{\link{kfs}} method. Otherwise, the model is returned as
	  is.}
	\item{phi}{a vector of initial values of hyperparamters. Note
	  that these have to be in the right order. Best advice is to
	  leave this option to be \code{NULL} and then inspect the
	  returned result using \code{phi(ssm)}. }
	\item{m0}{a vector with the initial state vector.}
	\item{C0}{a matrix with the variance matrix of the initial state.}
	\item{Fmat}{a function giving the regression matrix at a given
	  timepoint.}
	\item{Gmat}{a function giving the evolution matrix at a given
	  timepoint.} 
	\item{Wmat}{a function giving the evolution variance matrix at a given
	  timepoint.} 
	\item{Vmat}{a function giving the observation variance matrix at
	  a given timepoint.}
	\item{value}{an object to be assigned to the element of the
	  state space model.}
}
%\details{
%}
\value{
  An object of class \code{ssm} with the following components
  \item{ss}{an object of class \code{\link{SS}} describing the state
    space model. In addition, the \code{ss} object contains the
    components \code{family} and \code{ntotal} (for binomial case).}
}
\author{Claus Dethlefsen and \enc{Søren}{Soren} Lundbye-Christensen.}
\seealso{\code{\link{SS}}, \code{\link{extended}}}
\examples{
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
}
\keyword{models}% at least one, from doc/KEYWORDS
