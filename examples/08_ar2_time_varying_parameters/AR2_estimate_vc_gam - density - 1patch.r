library(cyclar)

# SIMULATION
# Evaluate linear change in parameters:
n <- 500
time <- ts(1:n)

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
#a0 <- 0
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

################################################################################
# SIMULATION

a.tv <- function(x,a,b,c) a / (1 + exp(-(x - b) / -c))

curve(a.tv(x,a=1,b=-10,c=100), from=0, to=500)

curve(a.tv(x,a=-1,b=-10,c=100), from=0, to=500)

n <- 500
a1.t <- a.tv(x=1:n, a=1, b=-10, c=100) # Direct d-d
a2.t <- a.tv(x=1:n, a=-1, b=-10, c=100/2) # Delayed d-d
a0 <- 2.1
a0.t <- rep(a0,n)

# Plot changes in parameter values over time:
par(mfrow=c(1,2))
plot(1:n, a1.t, type="l", col=2, lwd=2, xlab="Time", ylab="a2", font.lab=2, las=1,
main="Time-varying parameter: a1", bty="l")
plot(1:n, a2.t, type="l", col=2, lwd=2, xlab="Time", ylab="a2", font.lab=2, las=1,
main="Time-varying parameter: a2", bty="l")
par(mfrow=c(1,1))

ar2.plot(xlim=c(-2,2),ylim=c(-1,0),triangle=FALSE)
ar2.arrows(a1=a1.t, a2=a2.t)

a1.t <- a1.t -1

################################################################################
# DATA GENERATION & PARAMETER ESTIMATION

npatch <- 100
z.const <- sapply(1:npatch, function(i) ar2.sim(n=n, mu=a0.t[1], a1=a1.t[1], a2=a2.t[1], sd=sqrt(0.2), plot=TRUE))
z <- sapply(1:npatch, function(i) ar2.sim.tv(n=n, a0=a0.t, a1=a1.t, a2=a2.t, sd=sqrt(0.2), plot=TRUE))

plot(x=c(-2,6), y=c(0,.75), type="n",xlab="Starting value", ylab="Density", font.lab=2, las=1, main="Start values for simulations")
lines(density(z[1,]),lwd=2)
abline(v=mean(z[1,]),col=1,lty=2)
lines(density(z.const[1,]),col=2,lwd=2)
abline(v=mean(z.const[1,]),col=2,lty=2)
legend("topright", c("VCM","Constant"), col=c(1,2), lwd=c(2,2), bty="n", cex=1.2)

# Plot constant model
plot(x=c(0,n),y=range(z.const), type="n", xlab="Time", ylab="log(Population density)", font.lab=2, las=1,bty="l",
	main="Simulated time series, Constant")
abline(h=a0,lty=2)
for (i in 1:ncol(z)) points(z.const[,i], type="l", col=i)

# Plot VCM model
plot(x=c(0,n),y=range(z), type="n", xlab="Time", ylab="log(Population density)", font.lab=2, las=1,bty="l",
	main="Simulated time series, VCM")
abline(h=a0,lty=2)
for (i in 1:ncol(z)) points(z[,i], type="l", col=i)

################################################################################
# Construct data for parameter estimation:

dat <- lapply(1:ncol(z), function(i) growth.rate(exp(z[,i]),lags=2))
dat <- lapply(1:ncol(z), function(i) na.omit(dat[[i]]))
sapply(1:ncol(z), function(i) nrow(dat[[i]])) # Check number of observations per patch

i <- sample(size=1, npatch, replace=FALSE) # Select 1 of npatch simulated datasets!
dat <- growth.rate(exp(z[,i]),lags=2)
dat <- na.omit(dat)

time <- start(dat)[1]:end(dat)[1]
b.gam <- gam(Lag0 ~ s(time,bs="ps",m=c(2,2),by=Lag1) + s(time,bs="ps",m=c(2,2),by=Lag2) + s(time), data=dat)

summary(b.gam)
par(mfrow=c(1,3))
	plot(b.gam)
par(mfrow=c(1,1))

gam.check(b.gam)

b.gam.fit <- ar2.vc.gam.fit(b.gam,method="raw")

ar2.plot()
ar2.arrows(a1.t, a2.t, lwd=1)
ar2.arrows(b.gam.fit[[1]]$fit, b.gam.fit[[2]]$fit, lwd=1, col=4)

par(mfrow=c(1,2))

plot (x=c(0,n), y=c(-2,2), type="n", xlab="Time", ylab="a1", main="Time-varying parameter: a1", bty="l", font.lab=2, las=1)
points(time, b.gam.fit[[1]][,"fit"],col=1,pch=1,type="l")
points(time, a1.t[b.gam.fit[[1]][,"x"]],col=2,type="l")

plot (x=c(0,n), y=c(-2,2), type="n", xlab="Time", ylab="a2", main="Time-varying parameter: a2", bty="l", font.lab=2, las=1)
points(time, b.gam.fit[[2]][,"fit"],col=1,pch=1,type="l")
points(time, a2.t[b.gam.fit[[2]][,"x"]],col=2,type="l")

par(mfrow=c(1,1))
