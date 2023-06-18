library(cyclar)

# CAUTION! Code does not work - not yet updated.

p <- as.numeric(ar2.parms(k=4, v=2))
p <- as.numeric(ar2.parms(k=7, v=2))
# p <- c(0.7,0.1)

ar2.dyn(p)

N <- ar2.sim(n=50, a0=0, a1=p[1], a2=p[2], sd=sqrt(0.2))
gr.N <- growth.rate(N$exp.ts.x, log=TRUE)

par(mfrow=c(1,2))
plot(N$ts.x,xlab="Time",ylab="Population density",font.lab=2,las=1,main="Population density")
plot(gr.N[,"R"],xlab="Time",ylab="R",font.lab=2,las=1,main="Population growth rate")
abline(h=0, lty=2)
par(mfrow=c(1,1))

lag.plot(N$ts.x,lag=2,do.lines=FALSE, main="Population density")

lag.plot(na.omit(gr.N[,"R"]),lag=2,do.lines=FALSE, main="Population growth rate")

fma <- lm(as.numeric(gr.N[,"R"]) ~ as.numeric(gr.N[,"Lag0"]))
fmb <- lm(as.numeric(gr.N[,"R"]) ~ as.numeric(gr.N[,"Lag1"]))

par(mfrow=c(2,1))
plot(as.numeric(gr.N[,"Lag0"]), as.numeric(gr.N[,"R"]), xlab=expression(x[t]), ylab="R", las=1)
abline(fma)
plot(as.numeric(gr.N[,"Lag1"]), as.numeric(gr.N[,"R"]), xlab=expression(x[t-1]), ylab="R", las=1)
abline(fmb)
par(mfrow=c(1,1))

#plot(resid(fma), na.omit(as.numeric(gr.N[,"R"])), xlab=expression(x[t-1]), ylab="R", las=1)

par(mfrow=c(1,2))
acf(N$ts.x, main="Population density")
acf(gr.N[,"R"], na.action=na.omit, main="Population growth rate")
par(mfrow=c(1,1))

acf.N <- N$acf$acf
acf.gr.N <- acf(na.omit(gr.N[,"R"]),plot=FALSE)$acf

#plot(acf.N, acf.gr.N)

acf2AR(acf.N)[2,1:2]
acf2AR(acf.gr.N)[2,1:2]
p


growth.rate(N$ts.x, log=FALSE)[,"R"]
growth.rate(N$exp.ts.x, log=TRUE)[,"R"]
diff(N$ts.x)

ar.univar(x=N$ts.x, log=FALSE, constant=0, aic=FALSE, order.max=2, method="arima")$ar
ar.univar(x=gr.N[,"R"], log=FALSE, constant=0, aic=FALSE, order.max=2, method="arima")$ar

ar.univar(x=N$ts.x, log=FALSE, constant=0, aic=FALSE, order.max=2, method="ar.ols")$ar

###
# Centering variables (give the same results...)
lag.mat <- lags(N$ts.x)
lag.mat[,"Lag0"]<-scale(lag.mat[,"Lag0"])
lag.mat[,"Lag1"]<-scale(lag.mat[,"Lag1"])
lag.mat[,"Lag2"]<-scale(lag.mat[,"Lag2"])
apply(lag.mat,2,sd)
apply(lag.mat,2,mean)

ar2.gls(x=lag.mat, method="raw")
ar2.gls(x=lag.mat, method="raw")$coef

ar.univar(x=N$ts.x[3:50],log=FALSE)$coef
ar.univar(x=lag.mat[,"Lag0"],log=FALSE)$coef
###

ar2.lm(x=gr.N, method="raw")$coefs
ar2.lm(x=gr.N, method="growth rate")$coefs

ar2.lm(x=gr.N, method="raw")$ictab
ar2.lm(x=gr.N, method="growth rate")$ictab

fm1 <- ar2.lm(x=gr.N, method="raw")
fm2 <- ar2.lm(x=gr.N, method="growth rate")

fm1$mod
fm1$coefs

summary(fm1$mod$fm12)
summary(fm2$mod$fm12)

confint(fm1$mod$fm12)
confint(fm2$mod$fm12)

summary(fm1$mod$fm12)$r.squared
summary(fm2$mod$fm12)$r.squared

ar.univar(x=N$ts.x, log=FALSE, constant=0, aic=FALSE, order.max=2, method="ar.ols")$ar
summary(fm1$mod$fm12)

acf(residuals(fm1$models$fm12))

# How is sigma^2 estimated in ar.ols???
ar2.plot()
# points(ar2.coef.fox[,1], ar2.coef.fox[,2], pch=rep(c(16,15),each=c(6,3)), col=c(rep(2,6),rep(4,3)), cex=1.3)
pcoef1 <- ar2.lm(x=gr.N, method="raw")$coefs[4,]
pcoef2 <- ar2.lm(x=gr.N, method="growth rate")$coefs[4,]

pcoef3 <- ar.univar(x=N$ts.x, log=FALSE, constant=0, aic=FALSE, order.max=2, method="ar.ols")$coef
pcoef4 <- ar.univar(x=na.omit(gr.N[,"R"]), log=FALSE, constant=0, aic=FALSE, order.max=2, method="ar.ols")$coef

points(p[1], p[2], pch=1, col=3, cex=1.3)
points(pcoef1[2], pcoef1[3], pch=15, col=4, cex=1.3)
points(pcoef2[2]+1, pcoef2[3], pch=16, col=2, cex=1.3)

points(pcoef3[1,1], pcoef3[1,2], pch=16, col=5, cex=1.3)
points(pcoef4[1,1], pcoef4[1,2], pch=16, col=6, cex=1.3)


x <- arima.sim(list(order=c(2,0,0), ar = c(a1,a2)), n=n, sd=sqrt(ar.var), n.start=burnin, rand.gen=rnorm) + a0

fm.r.lm <- ar2.lm(growth.rate(exp(x),log=TRUE),method="growth rate")
summary(fm.r.lm$models$fm12)
fm.x.lm <- ar2.lm(growth.rate(exp(x),log=TRUE),method="raw")
summary(fm.x.lm$models$fm12)

fm.r.glm <- ar2.glm(growth.rate(exp(x),log=TRUE),method="growth rate")
summary(fm.r.glm$models$fm12)
fm.x.glm <- ar2.glm(growth.rate(exp(x),log=TRUE),method="raw")
summary(fm.x.glm$models$fm12)

fm.r.lm$ictab
fm.x.lm$ictab

fm.r.glm$ictab
fm.x.glm$ictab
