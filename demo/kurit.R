data(kurit)  ## West & Harrison, page 40

m1 <- SS(y=kurit,
         Fmat=function(tt,x,phi) return(matrix(1)),
         Gmat=function(tt,x,phi) return(matrix(1)),
         Wmat=function(tt,x,phi) return(matrix(5)),
         Vmat=function(tt,x,phi) return(matrix(100)),
         m0=matrix(130),C0=matrix(400)
         )

plot(m1$y)
m1.f <- kfilter(m1)
m1.s <- smoother(m1.f)
lines(m1.f$m,lty=2,col=2)
lines(m1.s$m,lty=2,col=2)

m2 <- m1

Wmat(m2) <- function(tt,x,phi) {
  if (tt==10) return(matrix(900))
  else return(matrix(5))
}

m2.f <- kfilter(m2)
m2.s <- smoother(m2.f)
lines(m2.f$m,lty=2,col=4)
lines(m2.s$m,lty=2,col=4)
