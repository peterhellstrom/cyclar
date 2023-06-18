library(cyclar)

# SIMULATION
# Evaluate linear change in parameters:
n <- 50
time <- ts(1:n)

##
a1.intercept <- 0.7
a1.slope <- -0.03

a2.intercept <- -0.5
a2.slope <- -0.0075

##
a1.intercept <- -.8
a1.slope <- 0.035

a2.intercept <- -0.8
a2.slope <- +0.0075

##
a1.t <- a1.intercept + (a1.slope*time)
a2.t <- a2.intercept + (a2.slope*time)

a0 <- 2.1
a0 <- 0
a0.t <- rep(a0,(n))

# Plot changes in parameter values over time:
par(mfrow=c(1,2))

plot (x=c(0,n), y=c(-2,2), type="n", xlab="Time", ylab="a1", main="Time-varying parameter: a1", bty="l", font.lab=2, las=1)
abline(a=a1.intercept, b=a1.slope, col=2, lwd=2)

plot (x=c(0,n), y=c(-2,2), type="n", xlab="Time", ylab="a2", main="Time-varying parameter: a2", bty="l", font.lab=2, las=1)
abline(a=a2.intercept, b=a2.slope, col=2, lwd=2)
par(mfrow=c(1,1))

ar2.plot()
ar2.arrows(a1.t,a2.t)

# SIMULATION

a.tv <- function(x,a,b,c) a / (1 + exp(-(x - b) / -c))

n <- 50
a1.t <- a.tv(x=1:n, a=1.3, b=-0.0025, c=10) # Direct d-d
a2.t <- a.tv(x=1:n, a=-0.8, b=30, c=8) # Delayed d-d
a0 <- 2.1
a0 <- 0
a0.t <- rep(a0,n)

# Plot changes in parameter values over time:
par(mfrow=c(1,2))
plot(1:n, a1.t, type="l", col=2, lwd=2, xlab="Time", ylab="a2", font.lab=2, las=1,
main="Time-varying parameter: a1", bty="l")
plot(1:n, a2.t, type="l", col=2, lwd=2, xlab="Time", ylab="a2", font.lab=2, las=1,
main="Time-varying parameter: a2", bty="l")
par(mfrow=c(1,1))

ar2.plot(xlims=c(-2,2),ylims=c(-1,0),triangle=FALSE)
ar2.arrows(a1=a1.t, a2=a2.t)

# DATA GENERATION & PARAMETER ESTIMATION

npatch <- 100
z <- sapply(1:npatch, function(i) ar2.sim.tv(n=50, a0=a0.t, a1=a1.t, a2=a2.t, sd=sqrt(0.2))$ts.x)
z <- ts(z, start=1)

# Plot
plot(x=c(0,n),y=range(z), type="n", xlab="Time", ylab="log(Population density)", font.lab=2, las=1,bty="l",main="Simulated time series")
abline(h=a0,lty=2)
for (i in 1:ncol(z)) points(z[,i], type="l", col=i)

# Plot changes in variability
z.s <- swin.var(x=z, swin=6, stat=sd)
matplot(z.s,type="l",xlab="Time (center of time-window)",ylab="S-index",main="Changes in crude variability over time", font.lab=2, las=1, bty="l")

# Construct data for parameter estimation:

dat <- lapply(1:ncol(z), function(i) growth.rate(exp(z[,i]),lags=2))
dat <- lapply(1:ncol(z), function(i) na.omit(dat[[i]]))
sapply(1:ncol(z), function(i) nrow(dat[[i]])) # Check number of observations per patch

vc.fit <- lapply(1:npatch, function(i) {
	time <- start(dat[[i]])[1]:end(dat[[i]])[1]
	b.gam <- gam(Lag0 ~ s(time,bs="ps",m=c(2,2),by=Lag1) + s(time,bs="ps",m=c(2,2),by=Lag2) + s(time), data=dat[[i]])
	b.gam.fit <- ar2.vc.gam.fit(b.gam,method="raw")
	b.gam.fit
})

ar.fit <- sapply(1:npatch, function(i) {
	ar.mod <- ar.univar(x=z[,i], log=FALSE, constant=0, method="arima", d=0, q=0, plot=FALSE)
	ar.mod.coef <- ar.mod$coef
	ar.mod.coef[1,]
})
ar.fit <- t(ar.fit)


ar2.plot()
for(i in 1:npatch) ar2.arrows(vc.fit[[i]][[1]][,"fit"], vc.fit[[i]][[2]][,"fit"], col=i, length=0)
ar2.arrows(a1.t, a2.t, lwd=2)

ar2.plot()
for(i in 1:npatch) points(ar.fit[i,1],ar.fit[i,2],col=4,pch=16,cex=0.8)
ar2.arrows(a1.t, a2.t, lwd=1)

par(mfrow=c(1,2))

plot (x=c(0,n), y=c(-2,2), type="n", xlab="Time", ylab="a1", main="Time-varying parameter: a1", bty="l", font.lab=2, las=1)

for (i in 1:npatch) {
	time <- start(dat[[i]])[1]:end(dat[[i]])[1]
	b.gam.fit <- vc.fit[[i]]
	points(time, b.gam.fit[[1]][,"fit"],col=i,pch=1,type="l")
}
points(time, a1.t[b.gam.fit[[1]][,"x"]],col=2,pch=16,cex=1.1)

plot (x=c(0,n), y=c(-2,2), type="n", xlab="Time", ylab="a2", main="Time-varying parameter: a2", bty="l", font.lab=2, las=1)
for (i in 1:npatch) {
	time <- start(dat[[i]])[1]:end(dat[[i]])[1]
	b.gam.fit <- vc.fit[[i]]
	points(time, b.gam.fit[[2]][,"fit"],col=i,pch=1,type="l")
}
points(time, a2.t[b.gam.fit[[2]][,"x"]],col=2,pch=16,cex=1.1)

par(mfrow=c(1,1))

