################################################################################
library(cyclar)
################################################################################

p <- ar2.parms(k=5,v=2)
ar2.period(p$a1,p$a2)

n <- 100
nrep <- 200
Ps <- matrix(nrow=nrep, ncol=4)

ar2.plot()

for (i in 1:nrep) {

	x <- ar2.sim(n=n, mu=0, a1=p$a1, a2=p$a2, sd=sqrt(0.2), plot=FALSE)

	# Compare model estimates:
	fm1 <- ar.univar(x=x, log=FALSE, constant=0, method="arima", d=0, q=0, plot=FALSE)$ar
	#fm4 <- ar.univar(x=x, log=FALSE, constant=0, method="arima", d=1, q=1, include.mean=TRUE)$ar
	#fm4 <- ar.univar(x=x, log=FALSE, constant=0, method="arima", d=1, q=1, include.mean=FALSE, plot=FALSE)$ar
	# Add first-order moving average term
	fm4 <- ar.univar(x=diff(x), log=FALSE, constant=0, method="arima", d=0, q=1, include.mean=FALSE, plot=FALSE)$ar

	p1 <- coef(fm1)
	p4 <- coef(fm4)

	Ps[i,1:2] <- p1[1:2]
	Ps[i,3:4] <- p4[1:2]

	points(p1[1], p1[2], col=2)
	points(p4[1], p4[2], col=4)
}


points(1+p$a1,p$a2,col=2,cex=1.5,pch=16)

par(mfrow=c(1,2))
	plot(Ps[,1],Ps[,3],main=expression(1 + beta[1]),xlab="AR",ylab="ARMA",font.lab=2,las=1)
		abline(0,1,lty=2,col=2)
	plot(Ps[,2],Ps[,4],main=expression(beta[2]),xlab="AR",ylab="ARMA",font.lab=2,las=1)
		abline(0,1,lty=2,col=2)
par(mfrow=c(1,1))


graphics.off()
