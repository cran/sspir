\name{forecast}
\alias{forecast}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Forecasts a Gaussian state space model}
\description{
  Forecasting of a Gaussian state space model.
}
\usage{
forecast(ss, k=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{an object of class \code{\link{SS}} defining the state space
  model which has been smoothed.}
  \item{k}{an integer specifying the number of time steps to be forecasted.}
}
\value{The object \code{ss} with updated components \code{m}, \code{C} and \code{mu}.}

\author{Anette Luther Christensen.}
\seealso{\code{\link{SS}},  \code{\link{smoother}}.}
\examples{
# Formulates Gaussian state space model
m1 <- SS(
     Fmat = function(tt,x,phi){matrix(cos(2*pi*tt/365))},
     Gmat = function(tt,x,phi){matrix(1)},
     Wmat = function(tt,x,phi){matrix(phi)},
     Vmat = function(tt,x,phi){matrix(0.1)},
     phi  = c(1e-1),
     m0   = matrix(1),
     C0   = matrix(100)
)

# Simulates observations
set.seed(984375)
m1 <- recursion(m1, 365)
plot(m1)

# Fits model
fit <- kfs(m1)
plot(fit)

# Change format of variances as these are time-invariant
fit$Wmat <- matrix(m1$phi)
fit$Vmat <- matrix(0.1)

# Estimates variances by use of EM algorithm
est.fit <- EMalgo(fit)


# Forecasting
fcast <- forecast(fit, k=100)

plot(fcast)
}
\keyword{models}% at least one, from doc/KEYWORDS
