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
\author{Claus Dethlefsen and Søren Lundbye-Christensen.}
\seealso{\code{\link{ssm}}, \code{\link{kfilter}},
  \code{\link{smoother}}, \code{\link{getFamily}}.}
\keyword{models}% at least one, from doc/KEYWORDS
