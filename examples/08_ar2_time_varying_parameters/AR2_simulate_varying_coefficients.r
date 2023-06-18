library(cyclar)

# Make parameters a function of time:
# Example function

a.tv <- function(x,a,b,c) a / (1 + exp(-(x - b) / -c))

n <- 50
a1.t <- a.tv(x=0:n, a=1.3, b=-0.0025, c=10) # Direct d-d
a2.t <- a.tv(x=0:n, a=-0.8, b=30, c=8) # Delayed d-d
a0 <- 2.1
a0 <- 0
a0.t <- rep(2.1,(n+1))

# Investigate properties of function:
bi <- seq(-10,40,5)
ci <- seq(1,20,2)


plot(x=c(0,50),y=c(-1,1),type="n",xlab="Time",ylab="Parameter value", font.lab=2, las=1, main="b")
for (i in 1:length(bi)) curve(a.tv(x=x, a=-0.8, b=bi[i], c=8), from=0, to=50, n=101, col=i, lwd=2, add=T)
legend("topleft", legend=bi, col=1:length(bi), lwd=2, bty="n", title="b =")


plot(x=c(0,50),y=c(-1,1),type="n",xlab="Time",ylab="Parameter value", font.lab=2, las=1, main="c")
for (i in 1:length(ci)) curve(a.tv(x=x, a=-0.8, b=30, c=ci[i]), from=0, to=50, n=101, col=i, lwd=2, add=T)
legend("topleft", legend=ci, col=1:length(ci), lwd=2, bty="n", title="c =")

# Plot changes in parameter values over time:
par(mfrow=c(1,2))
plot(0:n, a1.t, type="l", col=2, lwd=2, xlab="Time", ylab="a2", font.lab=2, las=1,
main="Time-varying parameter: a1", bty="l")
plot(0:n, a2.t, type="l", col=2, lwd=2, xlab="Time", ylab="a2", font.lab=2, las=1,
main="Time-varying parameter: a2", bty="l")
par(mfrow=c(1,1))

ar2.plot(xlims=c(-2,2),ylims=c(-1,0),triangle=FALSE)
ar2.arrows(a1=a1.t, a2=a2.t)


# Estimate ar-coefficients from time-varying simulations:
ar2.plot(ylims=c(-1,0),triangle=FALSE,winw=7,winh=5)
nrep <- 100

for (i in 1:nrep) {

	x <- ar2.sim.tv(n=50, a0=a0.t, a1=a1.t, a2=a2.t, sd=sqrt(0.2))$ts.x
	#x <- growth.rate(exp(x), log=TRUE, lags=1)[,"R"]
	mod <- ar.univar(x, log=FALSE, constant=0)
	points(mod$coef[1,1], mod$coef[1,2], col=4, pch=16, cex=1)
}

ar2.arrows(a1=a1.t, a2=a2.t)

######
npatch <- 10
z <- sapply(1:npatch, function(i) ar2.sim.tv(n=51, a0=a0.t, a1=a1.t, a2=a2.t, sd=sqrt(0.2))$ts.x)

# Plot
plot(x=c(0,n),y=range(z), type="n", xlab="Time", ylab="log(Population density)", font.lab=2, las=1,bty="l",main="Simulated time series")
abline(h=a0.t,lty=2)
for (i in 1:ncol(z)) points(z[,i], type="l", col=i)

plot(x=c(0,n),y=range(exp(z)), type="n", xlab="Time", ylab="Population density", font.lab=2, las=1,bty="l",main="Simulated time series")
abline(h=a0.t,lty=2)
for (i in 1:ncol(z)) points(exp(z[,i]), type="l", col=i)

# Plot changes in variability
z.s <- sapply(1:npatch, function(i) swin.var(x=z[,i], swin=6))
matplot(z.s,type="l",xlab="Time (center of time-window)",ylab="S-index",main="Changes in crude variability over time", font.lab=2, las=1, bty="l")


# Save one time series as x and use for estimation:

i <- 10
x <- ts(z[,i])

xd <- detr.ts(x=start(x)[1]:end(x)[1], y=x, f=10)
xd <- xd$y.lws


# Construct data for parameter estimation:

dat <- growth.rate(exp(x),lags=1)
dat <- na.omit(dat)
time <- start(dat[,"R"])[1]:end(dat[,"R"])[1]


dat <- data.frame(
R = as.numeric(dat[,"R"]),
Lag0 = as.numeric(dat[,"Lag0"]),
Lag1 = as.numeric(dat[,"Lag1"])
)

datd <- growth.rate(exp(xd),lags=1)
datd <- na.omit(datd)
time <- start(datd[,"R"])[1]:end(datd[,"R"])[1]


datd <- data.frame(
R = as.numeric(datd[,"R"]),
Lag0 = as.numeric(datd[,"Lag0"]),
Lag1 = as.numeric(datd[,"Lag1"])
)

################################################################################

	b.gam <- gam(R ~
			s(time, bs="ps", m=c(2,2), by=Lag0) +
			s(time, bs="ps", m=c(2,2), by=Lag1),
			data=dat)

	b.gamd <- gam(R ~
			s(time, bs="ps", m=c(2,2), by=Lag0) +
			s(time, bs="ps", m=c(2,2), by=Lag1),
			data=datd)

	b.gam.fit <- ar2.vc.gam.fit(b.gam, method="raw")
	bd.gam.fit <- ar2.vc.gam.fit(b.gamd, method="raw")

	ar2.plot()
	ar2.arrows(a1.t, a2.t, col=2, length=0.1)
	ar2.arrows(1 + b.gam.fit[[1]]$fit, b.gam.fit[[2]]$fit, col=4, length=0.1)
	ar2.arrows(1 + bd.gam.fit[[1]]$fit, bd.gam.fit[[2]]$fit, col=5, length=0.1)
	legend("topright",legend=c("Simulated","Estimated","Estimated, detrended"), col=c(2,4,5), bty="n", lwd=c(1,1,1))

	ar2.vc.gam.plot(x=b.gam.fit,fit="R")
	ar2.vc.gam.plot(x=bd.gam.fit,fit="R")

	par(mfrow=c(1,2))
	plot(b.gam.fit[[1]]$x, 1 + b.gam.fit[[1]]$fit,ylim=range(c(a1.t,1 + b.gam.fit[[1]]$fit)),col=4,
		xlab="Time", ylab="Parameter value", main="Direct d-d", font.lab=2, las=1)
	points(time,a1.t[time],col=2)
	points(bd.gam.fit[[1]]$x, 1 + bd.gam.fit[[1]]$fit, col=5)
	legend("topright",legend=c("Simulated","Estimated","Estimated, detrended"), col=c(2,4,5), pch=c(1,1,1), bty="n", lwd=c(1,1,1))

	plot(b.gam.fit[[2]]$x, b.gam.fit[[2]]$fit,ylim=range(c(a2.t,b.gam.fit[[2]]$fit)),col=4,
		xlab="Time", ylab="Parameter value", main="Delayed d-d", font.lab=2, las=1)
	points(time,a2.t[time],col=2)
	points(bd.gam.fit[[2]]$x, bd.gam.fit[[2]]$fit, col=5)
	legend("topright",legend=c("Simulated","Estimated","Estimated, detrended"), col=c(2,4,5), pch=c(1,1,1), bty="n", lwd=c(1,1,1))
	par(mfrow=c(1,1))






