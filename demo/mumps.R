data(mumpsdat)
time <- 1:nrow(mumpsdat)
attach(mumpsdat)

m1 <- glm(mumps~poly(time,2)+ factor(mo),family=quasipoisson,data=mumpsdat)
cat("Dispersion:",summary(m1)$dispersion,"\n")

m2 <- glm(mumps~poly(time,2)+ polytrig(time,12),family=quasipoisson,data=mumpsdat)
cat("Dispersion:",summary(m2)$dispersion,"\n")

m3 <- ssm( mumps ~ -1 + tvar(polytime(time,1)) +
                  tvar(polytrig(time,12,1)),
                  family=poisson(link=log),time=time) 
m3$ss$phi["epsilon"] <- 0
m3$ss$phi["polytime(time, 1)time0"] <- 0
m3$ss$phi["polytime(time, 1)time1"] <- 0.0005
m3$ss$phi["polytrig(time,12,1)"] <- 0.0001
diag(m3$ss$C0) <- 1

m3.fit <- kfs(m3)

m4 <- glm(mumps~offset(m3.fit$mu[1,]),family=quasipoisson)
cat("Dispersion:",summary(m4)$dispersion,"\n")


apC <- unlist(lapply(m3.fit$C,function(M) diag(M)[1]))
year <- 1928 + (time-1)/12

#postscript("mumps.pdf",width=10,height=6,horizontal=FALSE)
#pdf("mumps.pdf",width=10,height=6)
par(mfrow=c(3,1),mar=c(0,5.1,0,2.1),oma=c(6,0,5,0))
plot(year,mumpsdat[,4],type='l',ylab='Number of Cases',xlab='',axes=F)
lines(year, exp(m3.fit$m[1,]),type='l',lwd=2)
axis(2)
box()
plot(year,12*atan2(m3.fit$m[4,],m3.fit$m[3,])/(2*pi),type='l',ylim=c(2.7,5.3),ylab='Peak',xlab='',lwd=2,axes=F)
abline(h=3,lty=3)
abline(h=4,lty=3)
abline(h=5,lty=3)
axis(2,at=c(3,4,5),labels=c("Apr","May","Jun"))
box()

plot(year,exp(2*sqrt(m3.fit$m[3,]^2 + m3.fit$m[4,]^2)),type='l',ylab='PT-ratio',xlab='Years',ylim=c(0,12),lwd=2,axes=F)
abline(h=0,lty=3)
abline(h=5,lty=3)
abline(h=10,lty=3)
axis(2)
axis(1,at=seq(1925,1975,5))
box()

#dev.off()

