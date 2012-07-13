\name{EMalgo}
\alias{EMalgo}

\title{Estimation of variance matrices in a Gaussian state space model}

\description{Estimates variance matrices of the observation- and latent process in a Gaussian state space model given as input, using the EM algorithm.}

\usage{
EMalgo(	ss,
        epsilon   = 1e-06,
        maxiter   = 50,
        Vstruc    = function(V) { return(V) },
        Wstruc    = function(W) { return(W) }
      )
}

\arguments{
  \item{ss}{an object of class \code{\link{SS}}.}
  \item{epsilon}{a (small) positive numeric giving the tolerance of the maximum relative differences of \code{Vmat} and \code{Wmat} between iterations.}
  \item{maxiter}{a positive integer giving the maximum number of iterations to run.}
  \item{Vstruc}{a function specifying the structure of the variance matrix of the observation model if such is available.}
   \item{Wstruc}{a function specifying the structure of the variance matrix of the state model if such is available.}
}

\details{
The EM algorithm is a two-step iterative maximisation algorithm used to estimate unknown and constant covariance matrices \eqn{V} and \eqn{W} in a Gaussian state space model. Let the state space model be given as in \code{\link{SS}} with constant covariance matrices.

The covariance matrix \eqn{V} is estimated by maximisation of \eqn{\log\left(p\left(Y_{t}|\theta_{t},V\right)\right)}{log(p(Y_t|\theta_t, V))} conditional on all information, denoted \eqn{D_n}. This is equivalent to minimise the expectation of \deqn{\log(V)+\frac{1}{n}\sum_{t=1}^{n} (Y_{t}-\mu_{t})^2 V^{-1}}{log(V) + n^-1 \Sigma (Y_t - \mu_t)^2V^-1} conditional on \eqn{D_n}.

This yields
\deqn{
\textrm{E}\left[\log(V)+\frac{1}{n}\sum_{t=1}^{n} (Y_{t}-\mu_{t})^2 V^{-1}|D_n\right]
}{
E[log(V) + n^-1 \Sigma (Y_t - \mu_t)^2 V^-1 | D_n]
}
\deqn{
\quad=\log(V) + \frac{1}{n}\sum_{t=1}^{n}\left(\textrm{trace}(V^{-1}F_{t}^\top \widetilde{C}_tF_{t})+(Y_{t}-F_{t}^\top \widetilde{m}_t )^2 V^{-1}\right)
}{
log(V) + n^-1 \Sigma [ trace( V^-1 F_t^T s.C F_t ) + (Y_t - F_t^T s.m_t)^2]
}
\deqn{
\quad=\log(V)+\textrm{trace}\bigg[\frac{1}{nV}\sum_{t=1}^{n}\bigg(F_{t}^\top \widetilde{C}_tF_{t}+ \left(Y_{t}-F_{t}^\top \widetilde{m}_t\right)^2 \bigg)\bigg].
}{
log(V) +  trace[ (nV)^-1 \Sigma ( F_t^T s.C F_t  + (Y_t - F_t^T s.m_t)^2)].
}

This is minimised for
\deqn{\hat{V}=\frac{1}{n}\sum_{t=1}^{n}\bigg(F_{t}^\top \widetilde{C}_tF_{t}+ \left(Y_{t}-F_{t}^\top \widetilde{m}_t\right)^2 \bigg).}{V^* = n^-1 \Sigma [ F_t^T s.C F_t  + (Y_t - F_t^T s.m_t)^2].}

Similarly \eqn{W} is estimated by minimising the expectation of
\deqn{\log|W|+\frac{1}{n}\sum_{t=1}^{n} (\theta_{t}-G_t\theta_{t-1})(\theta_{t}-G_t\theta_{t-1})^\top W^{-1}}{log|W| + n^-1 \Sigma (\theta_t - G_t\theta_t-1)(\theta_t - G_t\theta_t-1)^T W^-1} conditional on \eqn{D_n}.

This yields
\deqn{
\textrm{E}\left[\log|W|+\frac{1}{n}\sum_{t=1}^{n} (\theta_{t}-G_t\theta_{t-1})(\theta_{t}-G_t\theta_{t-1})^\top W^{-1}|D_n\right]
}{
E[log|W| + n^-1 \Sigma (\theta_t - G_t\theta_t-1)(\theta_t - G_t\theta_t-1)^T W^-1 | D_n]
}
\deqn{
\quad=\log|W| + \frac{1}{n}\sum_{t=1}^{n}\left(\textrm{trace}(W^{-1}L_{t})+(\widetilde{m}_{t}-G_{t} \widetilde{m}_{t-1} )(\widetilde{m}_{t}-G_{t} \widetilde{m}_{t-1} )^\top W^{-1}\right)
}{
log|W| + n^-1 \Sigma [ trace( W^-1 L_t ) + (s.m_t - G_t s.m_t-1)(s.m_t - G_t s.m_t-1)^T]
}
\deqn{
\quad=\log|W|+\textrm{trace}\bigg[\frac{1}{n}W^{-1}\sum_{t=1}^{n}\bigg(L_{t} + \left(\widetilde{m}_{t}-G_{t} \widetilde{m}_{t-1}\right)\left(\widetilde{m}_{t}-G_{t} \widetilde{m}_{t-1}\right)^\top \bigg)\bigg].
}{
log|W| +  trace[ (nW)^-1 \Sigma ( L_t  + (s.m_t - G_t s.m_t-1)(s.m_t - G_t s.m_t-1)^T)].
}

This is minimised for
\deqn{\hat{W}=\frac{1}{n}\sum_{t=1}^{n}\bigg(L_{t} + \left(\widetilde{m}_{t}-G_{t} \widetilde{m}_{t-1}\right)\left(\widetilde{m}_{t}-G_{t} \widetilde{m}_{t-1}\right)^\top \bigg),}{W^* = n^-1 \Sigma [ L_t  + (s.m_t - G_t s.m_t-1)(s.m_t - G_t s.m_t-1)^T],}

where \eqn{\widetilde{m}_t, \widetilde{C}_t}{s.m_t, s.C_t} denotes smoothed means and covariances and \eqn{L_t=\textrm{Var}[\theta_{t}-G_t\theta_{t-1}|D_n]}{L_t=Var[\theta_{t}-G_t\theta_{t-1}|D_n]}.

}

\value{
  \item{ss}{the output from \code{\link{smoother}}.}
  \item{Vmat.est}{the estimate of the observation variance matrix, which is provided if the input of \code{Vstruc} is of class \code{function}, otherwise as input of \code{Vmat}.}
  \item{Wmat.est}{the estimate of the observation variance matrix, which is provided if the input of \code{Wstruc} is of class \code{function}, otherwise as input of \code{Wmat}.}
  \item{loglik}{maximum value of log likelihood function for each iteration.}
  \item{iteration}{number of iterations upon convergence.}
}

\author{Anette Luther Christensen}

\seealso{\code{\link{kfilter}}, \code{\link{smoother}}, \code{\link{recursion}}.}

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
          Wmat = matrix(c(0.01,0,0,0,0.1,0,0,0,0.1), nrow=3),
          Vmat = matrix(1),
          m0   = matrix(c(0,0,0), nrow=1),
          C0   = matrix(c(1,0,0,0,0.001,0,0,0,0.001), nrow=3, ncol=3)
        )

## Simulates 100 observation from m1
m1 <- recursion(m1, 100)

## Specifies the correlation structure of W:
Wstruc <- function(W){	
            W[1,2:3] <- 0
            W[2:3,1] <- 0
            W[2,2]   <- (W[2,2]+W[3,3])/2
            W[3,3]   <- W[2,2]
            W[2,3]   <- 0
            W[3,2]   <- 0
            return(W)
          }

## Estimstes variances and covariances of W by use of the EM algorithm
estimates <- EMalgo(m1, Wstruc=Wstruc)
estimates$Wmat.est

## Plots estimated model
plot(estimates$model)
}
\keyword{models}
