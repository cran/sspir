data(mumps)
index <- 1:length(mumps)
#attach(mumpsdat)

m1 <- glm(mumps~poly(index,2)+ factor(cycle(mumps)),family=quasipoisson)
cat("Dispersion:",summary(m1)$dispersion,"\n")

m2 <- glm(mumps~poly(index,2)+ polytrig(index,12),family=quasipoisson)
cat("Dispersion:",summary(m2)$dispersion,"\n")

phi.start <- c(0,0,0.0005,0.0001)

m3 <- ssm( mumps ~ -1 + tvar(polytime(index,1)) +
                  tvar(polytrig(index,12,1)),
                  family=poisson(link=log),
          phi=phi.start,
          C0 = diag(4)
          ) 

m3.fit <- getFit(m3)

m4 <- glm(mumps~offset(m3.fit$mu),family=quasipoisson)
cat("Dispersion:",summary(m4)$dispersion,"\n")


apC <- unlist(lapply(m3.fit$C,function(M) diag(M)[1]))

#postscript("mumps.pdf",width=10,height=6,horizontal=FALSE)
#pdf("mumps.pdf",width=10,height=6)
par(mfrow=c(3,1),mar=c(0,5.1,0,2.1),oma=c(6,0,5,0))
plot(mumps,type='l',ylab='Number of Cases',xlab='',axes=FALSE)
lines(exp(m3.fit$m[,1]),type='l',lwd=2)
axis(2)
box()
plot(12*atan2(m3.fit$m[,4],m3.fit$m[,3])/(2*pi),type='l',ylim=c(2.7,5.3),ylab='Peak',xlab='',lwd=2,axes=FALSE)
abline(h=3,lty=3)
abline(h=4,lty=3)
abline(h=5,lty=3)
axis(2,at=c(3,4,5),labels=c("Apr","May","Jun"))
box()

plot(exp(2*sqrt(m3.fit$m[,3]^2 + m3.fit$m[,4]^2)),type='l',ylab='PT-ratio',xlab='Years',ylim=c(0,12),lwd=2,axes=FALSE)
abline(h=0,lty=3)
abline(h=5,lty=3)
abline(h=10,lty=3)
axis(2)
axis(1,at=seq(1925,1975,5))
box()

#dev.off()

