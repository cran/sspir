forecast <- function(ss, k=10){
ss$k <- k

m <- matrix(NA, nrow=ss$n+k, ncol=ss$p)
m[1:ss$n, 1:ss$p] <- ss$m
C <- ss$C
f <- matrix(NA, nrow=ss$n+k, ncol=ss$d)
f[1:ss$n, 1:ss$d] <- ss$f
Q <- ss$Q
mu <- matrix(NA, nrow=ss$n+k, ncol=ss$p)
mu[1:ss$n, 1:ss$p] <- ss$mu

m[ss$n+1,] <- ss$Gmat(ss$n+1,ss$x,ss$phi)%*%matrix(m[ss$n,],ncol=1)
if(class(ss$Wmat)=="function") {
C[[ss$n+1]] <- ss$Gmat(ss$n+1,ss$x,ss$phi)%*%C[[ss$n]]%*%t(ss$Gmat(ss$n+1,ss$x,ss$phi))+ss$Wmat(ss$n+1,ss$x,ss$phi)
}else{
C[[ss$n+1]] <- ss$Gmat(ss$n+1,ss$x,ss$phi)%*%C[[ss$n]]%*%t(ss$Gmat(ss$n+1,ss$x,ss$phi))+ss$Wmat
}
f[ss$n+1,] <- t(ss$Fmat(ss$n+1,ss$x,ss$phi))%*%matrix(m[ss$n+1,],ncol=1)
if(class(ss$Wmat)=="function") {
Q[[ss$n+1]] <- t(ss$Fmat(ss$n+1,ss$x,ss$phi))%*%C[[ss$n+1]]%*%ss$Fmat(ss$n+1,ss$x,ss$phi)+ss$Vmat(ss$n+1,ss$x,ss$phi)
}else{
Q[[ss$n+1]] <- t(ss$Fmat(ss$n+1,ss$x,ss$phi))%*%C[[ss$n+1]]%*%ss$Fmat(ss$n+1,ss$x,ss$phi)+ss$Vmat
}
mu[ss$n+1,] <- t(ss$Fmat(ss$n+1,ss$x,ss$phi))%*%matrix(m[ss$n,], ncol=1)


for(tt in 2:k){

m[ss$n+tt,] <- ss$Gmat(ss$n+tt,ss$x,ss$phi)%*%matrix(m[ss$n+tt-1,],ncol=1)
mu[ss$n+tt,] <- t(ss$Fmat(ss$n+tt,ss$x,ss$phi))%*%matrix(m[ss$n+tt,], ncol=1)
f[ss$n+1,] <- t(ss$Fmat(ss$n+tt,ss$x,ss$phi))%*%matrix(m[ss$n+tt-1,],ncol=1)
if(class(ss$Wmat)=="function") {
C[[ss$n+tt]] <- ss$Gmat(ss$n+tt,ss$x,ss$phi)%*%C[[ss$n+tt-1]]%*%t(ss$Gmat(ss$n+tt,ss$x,ss$phi))+ss$Wmat(ss$n+tt,ss$x,ss$phi)
} else {
C[[ss$n+tt]] <- ss$Gmat(ss$n+tt,ss$x,ss$phi)%*%C[[ss$n+tt-1]]%*%t(ss$Gmat(ss$n+tt,ss$x,ss$phi))+ss$Wmat
}
if(class(ss$Vmat)=="function") {
Q[[ss$n+tt]] <- t(ss$Fmat(ss$n+tt,ss$x,ss$phi))%*%C[[ss$n+tt]]%*%ss$Fmat(ss$n+tt,ss$x,ss$phi)+ss$Vmat(ss$n+tt,ss$x,ss$phi)
} else {
Q[[ss$n+tt]] <- t(ss$Fmat(ss$n+tt,ss$x,ss$phi))%*%C[[ss$n+tt]]%*%ss$Fmat(ss$n+tt,ss$x,ss$phi)+ss$Vmat
}


}

ss$m <- m
ss$C <- C
ss$mu <- mu
ss$f <- f
ss$Q <- Q

ss
}