\name{kfs}
\alias{kfs}
\alias{kfs.ssm}
\alias{kfs.SS}
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
  \code{\link{extended}} is called. The \code{kfs} function is called by
  default in the call to \code{\link{ssm} unless the option
    \code{fit=FALSE} is given.}
}
\value{The returned value from either \code{\link{smoother}} (Gaussian
  case) or \code{\link{extended}}.
}
\author{Claus Dethlefsen and \enc{Søren}{Soren} Lundbye-Christensen.}
\seealso{\code{\link{kfilter}}, \code{\link{smoother}},
  \code{\link{extended}}, \code{\link{ssm}}.}
\examples{
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
}
\keyword{models}% at least one, from doc/KEYWORDS
