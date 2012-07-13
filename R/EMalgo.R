EMalgo <- function(
	ss,
	epsilon = 1e-6,
	maxiter = 50,
	Vstruc  = function(V){ return(V) },
	Wstruc  = function(W){ return(W) }
	)
  {

  rdiffV <- 1000;rdiffW <- 1000;ite <- 0;loglike <- numeric();Vini <- ss$Vmat;Wini <- ss$Wmat;conv <- 200
  x <- ss$x
  while( conv > epsilon & ite <= maxiter ) {
	ite <- ite + 1
	filtered <- kfilter(ss)
	smoothed <- smoother(filtered)

	R <- matrix(nearPD(filtered$Gmat(1,x,phi) %*% filtered$C0 %*% t(filtered$Gmat(1,x,phi)) + filtered$Wmat)$mat, ncol=ss$p)
	B <- filtered$C0 %*% t(filtered$Gmat(1,x,phi)) %*% solve(R)
	smoothed$m0 <- t(t(filtered$m0) + B %*% (t(matrix(smoothed$m[1,],nrow=1)) - filtered$Gmat(1,x,phi) %*% t(filtered$m0)))
	smoothed$C0 <- matrix(nearPD(filtered$C0 + B %*% (smoothed$C[[1]] - R) %*% t(B))$mat, ncol=ss$p)

        step <- EMstep.sspir(filtered, smoothed, Vstruc = Vstruc, Wstruc = Wstruc)
	
	if(ite>1){
		conv <- (smoothed$loglik - loglike[ite-1])
	        }

	V <- step$V
	ss$Vmat <- V
	W <- step$W
	ss$Wmat <- W

	loglike[ite] <- smoothed$loglik
	}

Vest <- ss$Vmat
West <- ss$Wmat

ss$Vmat <- Vini;ss$Wmat <- Wini

list(model=ss,Vmat.est=Vest,Wmat.est=West,loglik=loglike,iterations=ite)
}

EMstep.sspir <- function( f,
                          s,
                          Vstruc = function(V){ return(V) },
                          Wstruc = function(W){ return(W) }
                        ){
  x <- f$x
  V <- t(s$Fmat(1,x,phi)) %*% s$C[[1]] %*% s$Fmat(1,x,phi) + (s$y[1]-t(s$Fmat(1,x,phi)) %*% t(matrix(s$m[1,],nrow=1))) %*% t(s$y[1]-t(s$Fmat(1,x,phi)) %*% t(matrix(s$m[1,],nrow=1)))

  R <- matrix(nearPD(s$Gmat(1,x,phi) %*% f$C0 %*% t(s$Gmat(1,x,phi)) + f$Wmat)$mat, ncol=f$p) #t
  B <- f$C0 %*% t(s$Gmat(1,x,phi)) %*% solve( R )             #t-1
  L <- matrix(nearPD(s$C[[1]]+ s$Gmat(1,x,phi) %*% s$C0 %*% t(s$Gmat(1,x,phi)) - s$C[[1]] %*% t( B ) %*% t(s$Gmat(1,x,phi)) - s$Gmat(1,x,phi) %*% B %*% t(s$C[[1]]))$mat, ncol=f$p)

  W <- matrix(nearPD(L+( (matrix(s$m[1,],ncol=1) - s$Gmat(1,x,phi) %*% t(s$m0)) %*% t(matrix(s$m[1,],ncol=1) - s$Gmat(1,x,phi) %*% t(s$m0) )) )$mat, ncol=f$p)

  for(t in 2:(s$n)){
    V <- V + t(s$Fmat(t,x,phi)) %*% s$C[[t]] %*% s$Fmat(t,x,phi) + (s$y[t]-t(s$Fmat(t,x,phi)) %*% t(matrix(s$m[t,],nrow=1))) %*% t(s$y[t]-t(s$Fmat(t,x,phi)) %*% t(matrix(s$m[t,],nrow=1)))

    R <- matrix(nearPD(s$Gmat(t,x,phi) %*% f$C[[t-1]] %*% t(s$Gmat(t,x,phi)) + f$Wmat)$mat, ncol=f$p) #t
    B <- f$C[[t-1]] %*% t(s$Gmat(t,x,phi)) %*% solve( R )             #t-1
    L <- matrix(nearPD(s$C[[t]]+ s$Gmat(t,x,phi) %*% s$C[[t-1]] %*% t(s$Gmat(t,x,phi)) - s$C[[t]] %*% t( B ) %*% t(s$Gmat(t,x,phi)) - s$Gmat(t,x,phi) %*% B %*% t(s$C[[t]]))$mat, ncol=f$p)
    W <- matrix(nearPD(L+( (matrix(s$m[t,],ncol=1) - s$Gmat(t,x,phi) %*% matrix(s$m[t-1,],ncol=1)) %*% t(matrix(s$m[t,],ncol=1) - s$Gmat(t,x,phi) %*% matrix(s$m[t-1,],ncol=1) )) )$mat, ncol=f$p)
  }

  V <- 1/s$n * V
  W <- 1/s$n * W

  V <- Vstruc(V)
  W <- matrix(nearPD(Wstruc(W))$mat, ncol=f$p)

  list(V=V,W=W)

}
