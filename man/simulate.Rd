\name{simulate}
\alias{simulate}
\alias{simulate.SS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Forwards filtering Backwards sampling}
\description{
  Draws \eqn{N} samples from the conditional distribution of
  \eqn{\theta_1,\ldots,\theta_n} given \eqn{y_1,\ldots,y_n} in a
  Gaussian state space model. 
}
\usage{
simulate(ss, N = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{an object of class \code{\link{SS}} with components \code{m}
    and \code{C} obtained from running \code{\link{kfilter}}.}
  \item{N}{an integer giving the number of simulated realizations wanted.}
}
\value{An array of dimension \eqn{p \times n \times N}. The \eqn{i}'th
  simulated state vector at time \eqn{t} is located in \code{[,t,i]}.
}
\references{
  C.K. Carter and R. Kohn (1994). On Gibbs Sampling for State Space
  Models. Biometrika, 3, 541-553.
}
\author{Claus Dethlefsen and Søren Lundbye-Christensen.}
\seealso{\code{\link{ssm}}, \code{\link{kfilter}},
  \code{\link{smoother}}.}
\keyword{models}% at least one, from doc/KEYWORDS
