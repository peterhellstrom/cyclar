library(cyclar)

# SIMULATION
npoints <- 8

coord.mat <- matrix(
c(
0,-0.2,
0.75, -0.2,
0.75, -0.8,
0, -0.8,
-1, -0.8,
-1, -0.2,
0.4, 0.3,
0,-0.5
),
ncol=2, nrow=npoints, byrow=TRUE, dimnames=list(1:npoints,c("x","y")))


route <- c(1,4,6,5,1)
# route <- c(4,8,1)
route <- c(1,2)
route <- c(2,8,5)
route <- c(5,1,3)

nsteps <- length(route)-1
yrs <- 100
#yrs <- 25

vals <- matrix(c(
	rep(yrs,nsteps),
	rep(0,nsteps),
	rep(0.2,nsteps)),
nrow=nsteps, ncol=3, byrow=FALSE)

ps <- ar2.sim.tv.gen(npoints=8, coord.mat=coord.mat, route=route, vals=vals, gam=TRUE)

a0.0 <- 2 # CHANGE THIS ONE, AFFECTS GAM FIT A LOT!!! MUCH WORSE FIT IF THIS VALUE >0!
a0 <- rep(a0.0, length(ps$a1))

################################################################################
# DATA GENERATION & PARAMETER ESTIMATION
n <- length(a0)
npatch <- 100
ar.var <- 0.2
z <- sapply(1:npatch, function(i) ar2.sim.tv(n=n, a0=a0, a1=ps$a1, a2=ps$a2, sd=sqrt(ar.var), plot.log=TRUE, vals=vals)$ts.x)
z <- ts(z, start=1)

# Plot
plot(x=c(0,n),y=range(z), type="n", xlab="Time", ylab="log(Population density)", font.lab=2, las=1,bty="l",main="Simulated time series")
abline(h=a0,lty=2)
for (i in 1:ncol(z)) points(z[,i], type="l", col=i)

# Plot changes in variability
z.s <- swin.var(x=z, swin=6, stat=sd, plot=TRUE)

matplot(z.s,type="l",xlab="Time (center of time-window)",ylab="S-index",main="Changes in crude variability over time", font.lab=2, las=1, bty="l")

################################################################################
# Construct data for parameter estimation:

dat <- lapply(1:ncol(z), function(i) growth.rate(exp(z[,i]),lags=2))
dat <- lapply(1:ncol(z), function(i) na.omit(dat[[i]]))
sapply(1:ncol(z), function(i) nrow(dat[[i]])) # Check number of observations per patch

################################################################################
# ADDED A SMOOTHER FOR TREND ALSO, without that one and a0.0 > 0, fit is totally WRONG!

vc.fit <- lapply(1:npatch, function(i) {

	time <- start(dat[[i]])[1]:end(dat[[i]])[1]

	b.gam <- gam(Lag0 ~
			s(time,bs="ps",m=c(2,2),by=Lag1) +
			s(time,bs="ps",m=c(2,2),by=Lag2) +
			s(time),
			data=dat[[i]])

	# summary(b.gam)
	# gam.check(b.gam)

	b.gam.fit <- ar2.vc.gam.fit(b.gam,method="raw")

	b.gam.fit
	})

ar.fit <- sapply(1:npatch, function(i) {
	ar.mod <- ar.univar(x=z[,i], log=FALSE, constant=0, method="arima", d=0, q=0, plot=FALSE)
	ar.mod.coef <- ar.mod$coef
	ar.mod.coef[1,]
	})
ar.fit <- t(ar.fit)

# Time-varying
ar2.plot()
for(i in 1:npatch) ar2.arrows(vc.fit[[i]][[1]][,"fit"], vc.fit[[i]][[2]][,"fit"], col=i, length=0.05)
ar2.arrows(ps$a1, ps$a2, lwd=2, length=0.3, end.only=TRUE)

# Time-invariant
ar2.plot()
for(i in 1:npatch) points(ar.fit[i,1],ar.fit[i,2],col=4,pch=16,cex=0.8)
ar2.arrows(ps$a1, ps$a2, end.only=TRUE, lwd=2, length=0.2)

# Plot each patch as a line over time
par(mfrow=c(1,2))

plot (x=c(0,n), y=c(-2,2), type="n", xlab="Time", ylab="a1", main="Time-varying parameter: a1", bty="l", font.lab=2, las=1)

for (i in 1:npatch) {
	time <- start(dat[[i]])[1]:end(dat[[i]])[1]
	b.gam.fit <- vc.fit[[i]]
	points(time, b.gam.fit[[1]][,"fit"],col=i,pch=1,type="l")
}
points(time, ps$a1[b.gam.fit[[1]][,"x"]],col=2,pch=16,cex=0.75)

plot (x=c(0,n), y=c(-2,2), type="n", xlab="Time", ylab="a2", main="Time-varying parameter: a2", bty="l", font.lab=2, las=1)
for (i in 1:npatch) {
	time <- start(dat[[i]])[1]:end(dat[[i]])[1]
	b.gam.fit <- vc.fit[[i]]
	points(time, b.gam.fit[[2]][,"fit"],col=i,pch=1,type="l")
}
points(time, ps$a2[b.gam.fit[[2]][,"x"]],col=2,pch=16,cex=0.75)

par(mfrow=c(1,1))

