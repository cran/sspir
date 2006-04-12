data(vandrivers)
#vd.time <- time(vandrivers$y)
vd <- ssm( y ~ tvar(1) + seatbelt + sumseason(time,12),
          family=poisson(link="log"),
          data=vandrivers,
          phi = c(1,0.0004945505),
          C0=diag(13)*1000
          )
vd <- getFit(vd)


mle <- function(phi.suggest,obj.ext) {
  cat(":",phi.suggest,":")
  obj.ext$phi[2] <- phi.suggest
  ext <- extended(obj.ext)
  lik <- ext$likelihood
  cat("loglik(",phi.suggest,") = ",lik)
  return( lik )
}

#res1 <- optim(par=4e-4,fn=mle,obj.ext=vd.res,method=c("L-BFGS-B"),
#              lower=0,upper=Inf,
#              control=list(parscale=1e-4,fnscale=-1),
#              hessian=FALSE)

#res1 <- optim(par=1,fn=mle,obj.ext=vd.res,method=c("SANN"),
#              control=list(fnscale=-1,maxiter=10),
#              hessian=FALSE)

#res1 <- optim(par=vd$ss$phi[2],fn=mle,method=c("Nelder-Mead"),
#              lower=rep(0,1),upper=rep(Inf,1),control=list(fnscale=-1),hessian=FALSE)
#res1 <- optim(par=vd$ss$phi[2],fn=mle,method=c("Nelder-Mead"),
#              control=list(trace=6,parscale=1e-4,fnscale=-1),hessian=FALSE)
#vd$ss$phi[2] <- res1$par


attach(vandrivers)

#pdf("vandrivers.pdf",width=10,height=6)
#postscript("vandrivers.pdf",width=10,height=6,horizontal=FALSE)
par(mfrow=c(1,1))
plot(y,ylim=c(0,20),ylab="Vandrivers killed",xlab="Time")
lines(exp(seatbelt*vd$m[,2] + vd$m[,1]),lwd=2)

vd.sd <- c()
for (i in 1:length(y)) {
  thisone <- vd$C[[i]][1:2,1:2]
  if (seatbelt[i]==0) { vd.sd <- c(vd.sd,sqrt(thisone[1,1])) }
  else
    vd.sd <- c(vd.sd,sqrt(sum(thisone)))
}

lines(exp(seatbelt*vd$m[,2] + vd$m[,1]+2*vd.sd),lty=2)
lines(exp(seatbelt*vd$m[,2] + vd$m[,1]-2*vd.sd),lty=2)
#dev.off()

cat("Reduction of casualties:",100*(1-exp(vd$m[1,2])),"%\n")
detach("vandrivers")
