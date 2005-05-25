\name{ssm}
\alias{ssm}
\alias{kfilter.ssm}
\alias{smoother.ssm}
\title{Define state-space model in a glm-style call.}
\description{
Use a glm-style formula and family arguments to setup a state space
model.
}
\usage{
ssm(formula, family = gaussian, data = list(), subset =
    NULL, time = NULL)
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
          a family function.  (See \code{\link{family}} and
	  \code{\link{getFamily}} for details of family functions.)} 
  \item{data}{an optional data frame containing the variables in the model.
           If not found in \code{data}, the variables are taken from
          \code{environment(formula)}, typically the environment from which
          \code{ssm} is called.}
  \item{subset}{an optional vector specifying a subset of observations to be
          used in the fitting process.}
  \item{time}{a vector giving the observation times, eg. \code{1:n}.}
}
%\details{
%}
\value{
  An object of class \code{ssm} with the following components
  \item{ss}{an object of class \code{\link{SS}} describing the state
    space model. In addition, the \code{ss} object contains the
    components \code{family} and \code{ntotal} (for binomial case).}
}
\author{Claus Dethlefsen and Søren Lundbye-Christensen.}
\seealso{\code{\link{SS}}, \code{\link{extended}}}
\examples{
data(vandrivers)
vd <- ssm( y ~ tvar(1) + seatbelt + sumseason(time,12),
          time=time, family=poisson(link="log"),
          data=vandrivers)
vd$ss$phi["(Intercept)"] <- exp(- 2*3.703307 )
vd$ss$C0 <- diag(13)*1000
vd.res <- kfs(vd)

plot( ts( t(vd.res$m[1:3,]) ))

attach(vandrivers)
plot(y,ylim=c(0,20))
lines(exp(vd.res$m[1,]+vd.res$m[2,]*seatbelt),lwd=2 )
}
\keyword{models}% at least one, from doc/KEYWORDS
