\name{kfs}
\alias{kfs}
\alias{kfs.ssm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{(Iterated extended) Kalman smoother}
\description{
  Runs the Kalman smoother on an \code{\link{ssm}} object if the
  family-attribute of the object is \code{gaussian}. If not, the
  iterated extended Kalman smoother is used.
}
\usage{
kfs(ss, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{an object of class \code{\link{ssm}}.}
  \item{\dots}{further arguments passed to the Kalman smoother.}
}
\details{
  The function is a wrapper for \code{\link{kfilter}} followed by
  \code{\link{smoother}} in the Gaussian case. For other family types,
  \code{\link{extended}} is called.
}
\value{The returned value from either \code{\link{smoother}} (Gaussian
  case) or \code{\link{extended}}.
}
\author{Claus Dethlefsen and Søren Lundbye-Christensen.}
\seealso{\code{\link{kfilter}}, \code{\link{smoother}},
  \code{\link{extended}}, \code{\link{ssm}}.}
\examples{
data(mumpsdat)
time <- 1:nrow(mumpsdat)
m3 <- ssm( mumps ~ -1 + tvar(polytime(time,1)) +
                  tvar(polytrig(time,12,1)),
                  family=poisson(link=log),time=time,data=mumpsdat) 
m3$ss$phi["epsilon"] <- 0
m3$ss$phi["polytime(time, 1)time0"] <- 0
m3$ss$phi["polytime(time, 1)time1"] <- 0.0005
m3$ss$phi["polytrig(time,12,1)"] <- 0.0001
diag(m3$ss$C0) <- 1

m3.fit <- kfs(m3)
}
\keyword{models}% at least one, from doc/KEYWORDS
