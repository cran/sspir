\name{recursion}
\alias{recursion}
\alias{recursion.SS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate from a Gaussian state space model}
\description{
  Draw a sample of \eqn{y_1,\ldots,y_n}, and
  \eqn{\theta_1,\ldots,\theta_n} from the state space model given as input.
}
\usage{
recursion(ss, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{an object of class \code{\link{SS}} defining the state space
  model.}
  \item{n}{an integer specifying the number of time steps to be simulated.}
}
\value{The object \code{ss} with updated components \code{y}, \code{mu}
  and \code{truetheta}.}
}
\author{Claus Dethlefsen and \enc{Søren}{Soren} Lundbye-Christensen.}
\seealso{\code{\link{ssm}}, \code{\link{kfilter}},
  \code{\link{smoother}}, \code{\link{getFamily}}.}
\examples{
data(kurit)
m1 <- SS(kurit)
phi(m1) <- c(100,5)
m0(m1) <- matrix(130)
C0(m1) <- matrix(400)

par(mfrow=c(2,1))
plot(recursion(m1,100))
phi(m1) <- c(5,100)
plot(recursion(m1,100))
}
\keyword{models}% at least one, from doc/KEYWORDS
