#pdf("gasdata.pdf",width=10,height=6)
#postscript("gasdata.pdf",width=10,height=6,horizontal=FALSE)
plot(log10(UKgas))
#dev.off()

phistart <- StructTS(log10(UKgas),type="BSM")$coef[c(4,1,2,3)]

#time <- 1:length(UKgas)
gasmodel <- ssm(log10(UKgas) ~ -1+tvar(polytime(time,1)) +
                tvar(sumseason(time,4)), 
                phi=phistart
                )
fit <- getFit(gasmodel)
colnames(fit$m) <- c("Trend","Slope","Season","Season2","Season3")

#pdf("gas.pdf",width=10,height=6)
#postscript("gas.pdf",width=10,height=6,horizontal=FALSE)

plot(fit$m[,1:3],main="")
abline(v=1971,lty=2)
abline(v=1979,lty=2)
#dev.off()

