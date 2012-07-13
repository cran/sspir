\name{Fkfilter}
\alias{Fkfilter}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kalman filtering of Gaussian state space model using\code{\link[KFAS]{KFS}}}
\description{
  Does Kalman filtering using \code{\link[KFAS]{KFS}}.
}
\usage{
Fkfilter(ss, tvar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{object of class \code{\link{SS}}.}
  \item{tvar}{vector with four elements indicating the dynamic nature of
  matrices in the object \code{SS}. The first element corresponds to
  \code{ss$Fmat}, the second to \code{ss$Gmat}, the third to
  \code{ss$Vmat} and the fourth to \code{ss$Wmat}. If \code{tvar[i]==1}
  the corresponding matrix is constant, otherwise \code{tvar[i]} is set
  to \code{ss$n}.}
}
\details{
  The Kalman filter implemented in KFAS provides faster computations than
  \code{\link{kfilter}}, however requires alternative specification of
  the state space model. This function transform an object of class
  \code{\link{SS}} into an object on which the filter implemented in KFAS
    may be applied.
}
\value{List with two objects. An object \code{kfas} with components as specified in \code{\link[KFAS]{KFS}}, as well as the input object \code{ss} with updated \code{m}, \code{C} and \code{loglik}.
}
\author{Anette Luther Christensen and Claus Dethlefsen.}
\seealso{\code{\link{Fkfs}}}
\examples{
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

}
\keyword{models}% at least one, from doc/KEYWORDS
