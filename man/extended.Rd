\name{extended}
\alias{extended}
\alias{extended.SS}
\alias{print.extended}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Iterated Extended Kalman Smoothing }
\description{
  An iterative procedure for calculation of the conditional mean and
  variance of the latent process in non-Gaussian state space models. The
  method calculates an approximating Gaussian state space model.
}
\usage{
extended( ss,
          maxiter = 50,
          epsilon = 1e-06,
          debug = FALSE
        )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{an object of class \code{\link{SS}}.}
  \item{maxiter}{a positive integer giving the maximum number of
    iterations to run.}
  \item{epsilon}{a (small) positive numeric giving the tolerance of the
    maximum relative differences of \eqn{m} and \eqn{C} between iterations.}
  \item{debug}{a logical. If \code{TRUE}, some extra information is printed.}
}

\details{This is the default method when using \code{\link{kfs}} on an
  object of class \code{\link{ssm}} when the \code{family} is not
  \code{gaussian}. The conditional mean and variance can be retrieved
  using \code{\link{getFit}} and are then stored in the attributes
  \code{m} and \code{C}, respectively.}

\references{
Durbin J, Koopman SJ (2001). Time series analysis by state space
methods. Oxford University Press.
}

\value{
  The object \code{ss} with updated components \code{m}, \code{C},
  \code{loglik}, \code{iteration}, \code{ytilde}, \code{x$vtilde},
  \code{mu}. These describe the approximating Gaussian state space model.
}
\author{Claus Dethlefsen and \enc{Søren}{Soren} Lundbye-Christensen.}
\seealso{\code{\link{ssm}}, \code{\link{kfilter}},
  \code{\link{smoother}}.}
\examples{
data(mumps)
index <- 1:length(mumps) # use 'index' instead of time
model <- ssm( mumps ~ -1 + tvar(polytime(index,1)),
              family=poisson(link=log))
results <- getFit(model)
plot(mumps,type='l',ylab='Number of Cases',xlab='',axes=FALSE)
lines( exp(results$m[,1]), lwd=2)
## Alternatives:
## results2 <- extended(model$ss)
## results3 <- kfs(model) ## yields the same
}
\keyword{models}% at least one, from doc/KEYWORDS
