\name{Specials}
\alias{polytrig}
\alias{polytime}
\alias{sumseason}
\alias{season}
\title{Special functions used in ssm formulas}
\description{
  Auxiliary function for modelling the time-dependent trend patterns
  when formulating \code{\link{ssm}} models.
}
\usage{
polytrig(time, period = 365.25, degree = 1)
polytime(time, degree = 1)
sumseason(time, period = 12)
season(time, period = 12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{a vector giving the observation times, eg. \code{1:n}. The
  variable 'time' is by default set to \code{time(y)}, where \code{y} is
the response in the formula in the \code{\link{ssm}} call. }
  \item{period}{the period of the seasonal pattern.}
  \item{degree}{the degree of the trigonometric polynomial or the
                degree of the polynomial time trend.}
	    }
\details{
\code{polytrig} defines a trigonometric polynomial in time,
\code{polytime} defines a polynomial time trend, \code{sumseason}
produces a seasonal pattern in which the factors sum to zero over the
period, and \code{season} defines a factor. The functions \code{polytrig} and
\code{season} can be used outside the \code{ssm} call. However,
\code{polytime} and \code{sumseason} are treated mainly inside the
\code{ssm} call.
}
\author{Claus Dethlefsen and \enc{Søren}{Soren} Lundbye-Christensen.}
\seealso{\code{\link{ssm}}, \code{\link{kfilter}},
  \code{\link{smoother}}, \code{\link{getFamily}}.}
\examples{
polytrig(1:10,degree=2)
season(1:12,period=4)
}
\keyword{models}% at least one, from doc/KEYWORDS
