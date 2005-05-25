#pdf("gasdata.pdf",width=10,height=6)
#postscript("gasdata.pdf",width=10,height=6,horizontal=FALSE)
plot(log10(UKgas))
#dev.off()

time <- 1:length(UKgas)
gasmodel <- ssm(log10(UKgas) ~ -1 + tvar(polytime(time,1)) +
                tvar(sumseason(time,4)), time=time)

gasmodel$ss$phi <- StructTS(log10(UKgas),type="BSM")$coef[c(4,1,2,3)]

fit <- kfs(gasmodel)

rownames(fit$m) <- c("Trend","Slope","Season","Season2","Season3")

#pdf("gas.pdf",width=10,height=6)
#postscript("gas.pdf",width=10,height=6,horizontal=FALSE)
plot( ts( t(fit$m[1:3,]),start=1960,frequency=4),main="" )
abline(v=1971,lty=2)
abline(v=1979,lty=2)
#dev.off()

