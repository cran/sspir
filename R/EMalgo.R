EMalgo <- function(
	ss,
	epsilon   = 1e-6,
	maxiter = 50,
	Vstruc = function(V){return(V)},
	Wstruc = function(W){return(W)},print.ite=F
	)
  {

  rdiffV <- 1000;rdiffW <- 1000;ite <- 1;loglike <- numeric();Vini <- ss$Vmat;Wini <- ss$Wmat
  x <- ss$x
  while( any(rdiffV > epsilon) & any(rdiffW > epsilon) & ite <= maxiter ) {
  if(print.ite==T)  print(ite)
	filtered <- kfilter(ss)
	smoothed <- smoother(filtered)

	R <- filtered$Gmat(1,x,phi) %*% filtered$C0 %*% t(filtered$Gmat(1,x,phi)) + filtered$Wmat
	B <- filtered$C0 %*% t(filtered$Gmat(1,x,phi)) %*% mysolve(R)
	smoothed$m0 <- t(t(filtered$m0) + B %*% (t(matrix(smoothed$m[1,],nrow=1)) - filtered$Gmat(1,x,phi) %*% t(filtered$m0)))
	smoothed$C0 <- filtered$C0 + B %*% (smoothed$C[[1]] - R) %*% t(B)
        
        step <- EMstep(filtered,smoothed,Vstruc=Vstruc,Wstruc=Wstruc)
	
	rdiffV <- abs(step$V - filtered$Vmat) %*% mysolve(filtered$Vmat)
	rdiffW <- abs(step$W - filtered$Wmat) %*% mysolve(filtered$Wmat)

	ss$Vmat <- step$V
	ss$Wmat <- step$W

	loglike[ite] <- smoothed$loglik
	ite <- ite + 1
	}

Vest <- ss$Vmat
West <- ss$Wmat

ss$Vmat <- Vini;ss$Wmat <- Wini

list(model=ss,Vmat.est=Vest,Wmat.est=West,loglik=loglike,iterations=ite)
}

EMstep <- function(f,s,Vstruc=function(V){ return(V)},Wstruc=function(W){ return(W)})
  {
x <- f$x
	V <- t(s$Fmat(1,x,phi)) %*% s$C[[1]] %*% s$Fmat(1,x,phi) +
	     (s$y[1]-t(s$Fmat(1,x,phi)) %*% t(matrix(s$m[1,],nrow=1))) %*%
	     t(s$y[1]-t(s$Fmat(1,x,phi)) %*% t(matrix(s$m[1,],nrow=1)))

  R <- s$Gmat(1,x,phi) %*% f$C0 %*% t(s$Gmat(1,x,phi)) + f$Wmat #t
	B <- f$C0 %*% t(s$Gmat(1,x,phi)) %*% mysolve( R )             #t-1
  L <- s$C[[1]]+ s$Gmat(1,x,phi) %*% s$C0 %*% t(s$Gmat(1,x,phi)) -
       s$C[[1]] %*% t( B ) %*% t(s$Gmat(1,x,phi)) -
       s$Gmat(1,x,phi) %*% B %*% t(s$C[[1]])

	W <- L+(t(matrix(s$m[1,],nrow=1)) - s$Gmat(1,x,phi) %*% t(s$m0)) %*% 
        t(t(matrix(s$m[1,],nrow=1)) - s$Gmat(1,x,phi) %*% t(s$m0))

    for(t in 2:(s$n)){

	V <- V + t(s$Fmat(t,x,phi)) %*% s$C[[t]] %*% s$Fmat(t,x,phi)+
	     (s$y[t]-t(s$Fmat(t,x,phi)) %*% t(matrix(s$m[t,],nrow=1))) %*%
	     t(s$y[t]-t(s$Fmat(t,x,phi)) %*% t(matrix(s$m[t,],nrow=1)))

  R <- s$Gmat(t,x,phi) %*% f$C[[t-1]] %*% t(s$Gmat(t,x,phi)) + f$Wmat #t
	B <- f$C[[t-1]] %*% t(s$Gmat(t,x,phi)) %*% mysolve( R )             #t-1
  L <- s$C[[t]]+ s$Gmat(t,x,phi) %*% s$C[[t-1]] %*% t(s$Gmat(t,x,phi)) -
       s$C[[t]] %*% t( B ) %*% t(s$Gmat(t,x,phi)) -
       s$Gmat(t,x,phi) %*% B %*% t(s$C[[t]])

	W <- W + L +
        (t(matrix(s$m[t,],nrow=1)) - s$Gmat(t,x,phi)%*%t(matrix(s$m[t-1,],nrow=1))) %*% 
       t(t(matrix(s$m[t,],nrow=1)) - s$Gmat(t,x,phi)%*% t(matrix(s$m[t-1,],nrow=1)))

	}

V <- 1/s$n * V
W <- 1/s$n * W


V <- Vstruc(V)
W <- Wstruc(W)

list(V=V,W=W)
	
}
