\name{kfilter}
\alias{kfilter}
\alias{kfilter.SS}
\alias{filterstep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kalman filter for Gaussian state space model}
\description{
  From an SS object, runs the Kalman filter to produce the conditional means and
  variances of the state vectors given the current time point.
}
\usage{
kfilter(ss)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{object of class \code{\link{SS}}.}
}
\details{
  The Kalman filter yields the distribution
  \deqn{(\theta_t|y_1,\ldots,y_t) \sim N(m_t, C_t)}{(\theta_t|y[,1:t]) ~ N(m_t, C_t)}
  through the recursion for \eqn{t=1,\ldots,n},
  \deqn{a_t = G_t  m_{t-1}}
  \deqn{R_t = G_t C_{t-1} G_t^T + W_t}
  \deqn{f_t = F_t^T  a_t}
  \deqn{Q_t = F_t^T  R_t  F_t + V_t}
  \deqn{e_t = y_t - f_t}
  \deqn{A_t = R_t  F_t  Q_t^{-1}}
  \deqn{m_t = a_t + A_t  e_t}
  \deqn{C_t = R_t - A_t  Q_t  A_t^T}
  Also, the log-likelihood is calculated.
}
\value{An object of class \code{\link{SS}} with the components \code{m},
  \code{C}, and \code{loglik} updated.
}
\author{Claus Dethlefsen and \enc{Søren}{Soren} Lundbye-Christensen.}
\seealso{\code{\link{SS}}, \code{\link{smoother}}}
\examples{
data(kurit)
m1 <- SS(kurit)
phi(m1) <- c(100,5)
m0(m1) <- matrix(130)
C0(m1) <- matrix(400)

m1.f <- kfilter(m1)
plot(m1$y)
lines(m1.f$m,lty=2)
}
\keyword{models}% at least one, from doc/KEYWORDS
