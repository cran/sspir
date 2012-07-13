
EMalgo.kfas <- function(
                        ss,
                        epsilon   = 1e-6,
                        maxiter   = 100,
                        Vstruc    = NULL,
                        Wstruc    = NULL,
                        print.ite = FALSE,
                        trace     = FALSE,
                        Wpath     = NULL,
                        Vpath     = NULL,
			T         = 2,
			S         = 2,
                        ...
                        )
{

    ite <- 0;loglike <- numeric();conv <- 1000

  m0 <- ss$m0;C0 <- ss$C0
  while( conv > epsilon & ite <= maxiter ) {
      ite <- ite + 1
      if(print.ite==TRUE)  print(ite)

  # Using KFAS to perform Kalman filtering and smoothing
	smoothed <- Fkfs(ss, ...)

    R <- smoothed$ss$Gmat(1, smoothed$ss$x, smoothed$ss$phi) %*% C0 %*% t(smoothed$ss$Gmat(1, smoothed$ss$x, smoothed$ss$phi)) + smoothed$ss$Wmat(1, smoothed$ss$x, smoothed$ss$phi)
	B <- C0 %*% t(smoothed$ss$Gmat(1, smoothed$ss$x, smoothed$ss$phi)) %*% solve(R)

	smoothed$ss$C0 <- C0 + B %*% (smoothed$ss$C[[1]] - R) %*% t(B)
    smoothed$ss$m0 <- t(t(m0) + B %*% (smoothed$ss$m[1,] - smoothed$ss$Gmat(1, smoothed$ss$x, smoothed$ss$phi) %*% t(m0)))

    step <- EMstep(smoothed$ss, m0, C0, T, S, if(is.null(Vstruc)){ V.est=no.strucV }else{ V.est=Vstruc }, if(is.null(Wstruc)){ W.est=no.strucW }else{ W.est=Wstruc })

  if(ite>1){
      conv <- (smoothed$ss$loglik - loglike[ite-1])
        }

      if(smoothed$ss$family$family=="gaussian"){V <- step$V}
      W <- step$W


	if(smoothed$ss$family$family=="gaussian"){ss$Vmat <- function(tt,x,phi) V}
	ss$Wmat <- function(tt,x,phi) W


        if(trace){
            write.table(matrix(step$W[lower.tri(step$W,diag=TRUE)], nrow=1), append = TRUE, col.names = FALSE, row.names = FALSE, file = Wpath)
            if(smoothed$ss$family$family=="gaussian"){
                write.table(matrix(step$V[lower.tri(step$V,diag=TRUE)], nrow=1), append = TRUE, col.names = FALSE, row.names = FALSE, file = Vpath)
            }
        }

	loglike[ite] <- smoothed$ss$loglik

	}

if(smoothed$ss$family$family=="gaussian"){Vest <- V}else{Vest <- NULL}
West <- W

  smoothed$ss$m0 <- m0;  smoothed$ss$C0 <- C0

  fit <- smoothed$ss
    if(conv<epsilon){cc <- TRUE}else{cc <- FALSE}

list(ss=smoothed$ss, Vmat.est=Vest, Wmat.est=West, loglik=loglike, iterations=ite, fit=fit, convergence=cc)
}

EMstep <- function(s, m0, C0, T, S, V.est=no.strucV, W.est=no.strucW){

    V.hat <- V.est(s)
    W.hat <- W.est(s, m0, C0, T, S)

    list(V=V.hat,W=W.hat)

}

no.strucV <- function(s){
    x <- s$x
    phi <- s$phi

    V <- t(s$Fmat(1,x,phi)) %*% s$C[[1]] %*% s$Fmat(1,x,phi) + (s$y[1] - t(s$Fmat(1,x,phi)) %*% s$m[1,]) %*% t(s$y[1] - t(s$Fmat(1,x,phi)) %*% s$m[1,])

    for(t in 2:(s$n)){
        	V <- V + t(s$Fmat(t,x,phi)) %*% s$C[[t]] %*% s$Fmat(t,x,phi) + (s$y[t]-t(s$Fmat(t,x,phi)) %*% t(matrix(s$m[t,],nrow=1))) %*% t(s$y[t] - t(s$Fmat(t,x,phi)) %*% t(matrix(s$m[t,],nrow=1)))
            }

    V <- 1/s$n*V

    return(V)
}

no.strucW <- function(s, m0, C0){
   x <- s$x
    phi <- s$phi

    R <- s$Gmat(1,x,phi) %*% C0 %*% t(s$Gmat(1,x,phi)) + s$Wmat(1,x,phi) #t
    B <- C0 %*% t(s$Gmat(1,x,phi)) %*% solve( R )             #t-1
    L <- s$C[[1]]+ s$Gmat(1,x,phi) %*% s$C0 %*% t(s$Gmat(1,x,phi)) -
       s$C[[1]] %*% t( B ) %*% t(s$Gmat(1,x,phi)) -
       s$Gmat(1,x,phi) %*% B %*% t(s$C[[1]])

    W1 <- L+(t(matrix(s$m[1,],nrow=1)) - s$Gmat(1,x,phi) %*% t(matrix(s$m0[1,],nrow=1))) %*%
        t(t(matrix(s$m[1,],nrow=1)) - s$Gmat(1,x,phi) %*% t(matrix(s$m0[1,],nrow=1)))


    for(t in 2:(s$n)){
        R <- s$Gmat(t,x,phi) %*% s$Cf[[t-1]] %*% t(s$Gmat(t,x,phi)) + s$Wmat(1,x,phi) #t
        B <- s$Cf[[t-1]] %*% t(s$Gmat(t,x,phi)) %*% solve( R )  #t-1
        L <- s$C[[t]]+ s$Gmat(t,x,phi) %*% s$C[[t-1]] %*% t(s$Gmat(t,x,phi)) - s$C[[t]] %*% t( B ) %*% t(s$Gmat(t,x,phi)) - s$Gmat(t,x,phi) %*% B %*% t(s$C[[t]])

        W1 <- W1 + L + (t(matrix(s$m[t,], nrow=1)) - s$Gmat(t,x,phi) %*% t(matrix(s$m[t-1,], nrow=1))) %*% t(t(matrix(s$m[t,], nrow=1)) - s$Gmat(t,x,phi) %*% t(matrix(s$m[t-1,], nrow=1)))

    }


    W <- 1/(2*s$n)*(W1 + t(W1))

    return(W)
}





