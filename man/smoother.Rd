\name{smoother}
\alias{smoother}
\alias{smoother.SS}
\alias{smootherstep}
\alias{smootherstep.uni}
\alias{print.Smoothed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Kalman smoother for Gaussian state space model}
\description{
  Based on the output from \code{\link{kfilter}}, this function runs the
  Kalman smoother to produce the conditional means and variances of the
  state vectors given all observations. 
}
\usage{
smoother(ss)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ss}{object of class \code{\link{SS}}, where the components
    \code{m} and \code{C} have been set by the \code{\link{kfilter}}.}
}
\details{
  The Kalman smoother yields the distribution
  \deqn{(\theta_t|y_1,\ldots,y_n) \sim N(m^*_t, C^*_t)}{(\theta_t|y[,1:n]) ~ N(m*_t, C*_t)}
  through the backward recursion for \eqn{t=n..1},
  \deqn{R_{t+1}= G_{t+1}  C_t  G_{t+1}^T + W_{t+1}}
  \deqn{B_t = C_t  G_{t+1}^T  R_{t+1}^{-1}}
  \deqn{m^*_t = m_t + B_t  ( m^*_{t+1} - G_{t+1}  m_t)}{m*_t = m_t + B_t  ( m*_{t+1} - G_{t+1}  m_t)}
  \deqn{C^*_t = C_t + B_t  ( C^*_{t+1} - R_{t+1} )B_t^T}{C*_t = C_t +
    B_t  ( C*_{t+1} - R_{t+1} ) B_t^T}
  where the matrices \eqn{F}, \eqn{G}, \eqn{V}, \eqn{W} are stored in
  the \code{SS} object as functions, eg. \code{Fmat(tt,x,phi)}, see
  \code{\link{SS}}. The vectors \eqn{m} and matrices \eqn{C} are set by
  the \code{\link{kfilter}} and are overwritten by the Kalman smoother.

  The smoother also calculates the signal, \eqn{\mu_t = F^T_t m^*_t}{\mu_t = F^T_t m*_t}.
}
\value{
  An object of class \code{\link{SS}} with the components \code{m},
  \code{C}, and \code{mu} updated.
}
\author{Claus Dethlefsen and Søren Lundbye-Christensen.}
\seealso{\code{\link{SS}}, \code{\link{kfilter}}}
\keyword{models}% at least one, from doc/KEYWORDS
