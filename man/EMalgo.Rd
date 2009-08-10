\name{EMalgo}
\alias{EMalgo}

\title{Estimation of variance matrices in a Gaussian state space model}

\description{Estimates variance matrices of the observation- and latent process in a Gaussian state space model given as input, using the EM algorithm.}

\usage{
EMalgo(ss, epsilon = 1e-06, maxiter = 50, Vstruc =
                 function(V) { return(V)}, Wstruc = function(W) {
                 return(W)}, print.ite = F)
}

\arguments{
  \item{ss}{an object of class \code{\link{SS}}.}
  \item{epsilon}{a (small) positive numeric giving the tolerance of the
    maximum relative differences of \code{Vmat} and \code{Wmat} between iterations.}
  \item{maxiter}{a positive integer giving the maximum number of
    iterations to run.}
	\item{Vstruc}{a function specifying the structure of the variance matrix of the observation model if such is available.}
	\item{Wstruc}{a function specifying the structure of the
	  variance matrix of the state model if such is available.}
	\item{print.ite}{A logical specifying whether information should
	  be printed after each iteration.}
}

\details{
The initial variance matrices are to be specified in the model specification and structures of the variance matrices may be specified by the user by the functions \code{Vstruc} and \code{Wstruc}. As default these are assigned \code{NA}, if not specified they are not estimated, hence it is possible to estimate only the observation variance matrix or the evolution variance matrix. The EM algorithm requires the variance matrices to be estimated are constant, however if a variance matrix is not be estimated it may be non-constant.

The output provided by the function is the smoothed Gaussian model along with the estimated variance matrices, maximum values of the log likelihood function for each iterations and the number of iteration upon convergence.
}

\value{
  \item{\code{ss}}{the value from \code{smoother}.}
  \item{\code{Vmat.est}}{the estimate of the observation variance matrix, which is provided if the input of \code{Vstruc} is of class \code{function}, otherwise as input of \code{Vmat}.}
  \item{\code{Wmat.est}}{the estimate of the observation variance matrix, which is provided if the input of \code{Wstruc} is of class \code{function}, otherwise as input of \code{Wmat}.}
  \item{\code{loglik}}{maximum value of log likelihood function for each iteration.}
  \item{\code{iteration}}{number of iterations upon convergence.}
}

\author{Anette Luther Christensen}

\seealso{\code{\link{kfilter}},\code{\link{smoother}},\code{\link{recursion}}.}

\examples{

m1 <- SS(

Fmat = function(tt,x,phi) { 
Fmat      <- matrix(NA,3,1)
Fmat[1,1] <- 1
Fmat[2,1] <- cos(2*pi*tt/12)
Fmat[3,1] <- sin(2*pi*tt/12) 
return(Fmat) },

Gmat = function(tt,x,phi) { 
return(matrix(c(1,0,0,0,1,0,0,0,1),nrow=3)) },

Wmat = matrix(c(0.01,0,0,0,0.1,0,0,0,0.1),nrow=3),
Vmat = matrix(1),
m0   = matrix(c(0,0,0),nrow=1),
C0   = matrix(c(1,0,0,0,0.001,0,0,0,0.001),nrow=3,ncol=3)
)

m1 <- recursion(m1,100)

Wstruc <- function(W){	
W[1,2:3]<-0
W[2:3,1]<-0
W[2,2]<-(W[2,2]+W[3,3])/2
W[3,3]<-W[2,2]
W[2,3]<-0
W[3,2]<-0
return(W)}

estimates <- EMalgo(m1,Wstruc=Wstruc)
plot(estimates$model)
}
\keyword{models}
