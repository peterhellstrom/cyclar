################################################################################
library(cyclar)
################################################################################

# Test the simulation function ar2.sim:

# Define parameters
k <- 5
v <- 2
p <- ar2.parms(k=k,v=v)

mu <- 2.1
a1 <- p$a1
a2 <- p$a2
ar2.var <- 0.2
n <- 200

# Note that var should be squared in order to get sd:
x <- ar2.sim(n=n, mu=mu, a1=a1, a2=a2, sd=sqrt(ar2.var), plot=TRUE)
ts.summary(x)
# expected sd
sqrt(ar2.ipv(a1,a2,sd=sqrt(ar2.var)))

################################################################################
# Simulate a "typical" North-Fennoscandian grey-sided vole pop:
# Periodicity 4.5 years, mean=2.1, "noise", var=0.2, ipv=2

k <- 4
v <- 3
mu <- 2.1 # Mean
p <- ar2.parms(k=k, v=v)
a1 <- p$a1
a2 <- p$a2
ar2.var <- 0.2

# Repeat simulations
# Number of replicate simulations
nreps <- 1000
ntime <- 1000

system.time(
sim.coef <- lapply(1:nreps, function(...) {

	# Simulate time series
	sim.ts <- ar2.sim(n=ntime, mu=mu, a1=a1, a2=a2, sd=sqrt(ar2.var), plot=FALSE)
	# Estimate AR2 coefficients from simulated time series and combine time series in a list
	yw.sim.ts <- ar.yw(sim.ts, order.max=2, aic=FALSE)

	list(
		ar = yw.sim.ts$ar,
		var.pred = yw.sim.ts$var.pred,
		x.mean = yw.sim.ts$x.mean,
		x = sim.ts)
	}
))

# Extract time series to a list
time.series <- lapply(1:nreps, function(i) sim.coef[[i]]$x)
# Extract all estimated coefficients in a list
coefs <- sapply(1:nreps, function(i) sim.coef[[i]]$ar)
ar2.dyn <- ar2.period(a1=coefs[1,]-1, a2=coefs[2,])

# Plot distribution of estimated coefficients
par(mfrow=c(1,2))
	hist(coefs[1,], breaks=30, freq=FALSE, las=1,
		main="Direct density-dependence", xlab=expression(paste("1+ ",beta,"1")), font.lab=2, col="lightgrey")
	abline(v=1+a1, col=2, lty=2)
	hist(coefs[2,], breaks=30, freq=FALSE, las=1,
		main="Delayed density-dependence", xlab=expression(paste(beta,"2")), font.lab=2, col="lightgrey")
	abline(v=a2, col=2, lty=2)
par(mfrow=c(1,1))

# Plot estimated coefficients in "triangle space"
ar2.plot(k=range(ar2.dyn[,"period"]), v=range(ar2.dyn[,"ipv"]), cex.contours=0.7)
points(coefs[1,], coefs[2,], col="red", cex=1, pch="+")
points(1+a1, a2, pch=21, col=1, bg="white", cex=1.5) # Input parameters for simulation

# Calculate average cycle length

# From fitted coefficients
ps <- t(sapply(1:nreps, function(i) sim.coef[[i]]$ar))
z <- ar2.period(a1=ps[,1]-1, a2=ps[,2], method="polyroot")
pairs(z)

par(mfrow=c(1,2))
hist(z[,"period"], breaks=30, freq=FALSE, col="lightgrey", xlab="Period", main="", font.lab=2, las=1)
abline(v=k, col=2, lty=2)
hist(z[,"ipv"], breaks=30, freq=FALSE, col="lightgrey", xlab="ipv", main="", font.lab=2, las=1)
abline(v=v, col=2, lty=2)
par(mfrow=c(1,1))

# From distance between peaks
period <- sapply(1:nreps, function(i) state.cycle(phase(sim.coef[[i]]$x,plot=FALSE))$mean.period)

hist(period, breaks=30, col="lightgrey", xlab="mean time between peaks", main="", font.lab=2, las=1)
abline(v=mean(period), lty=2)
abline(v=ar2.period(a1=a1,a2=a2)[,"period"], lty=2, col=2)

# Calculate s-index
s.index <- sapply(1:nreps, function(i) sd(time.series[[i]]))
hist(s.index, breaks=30, col="lightgrey", font.lab=2, las=1)

# Calculate max/min (ratio) for back-transformed densities
ratio <- sapply(1:nreps, function(x) {max(exp(time.series[[x]])) / min(exp(time.series[[x]]))})
hist(ratio, breaks=30, col="lightgrey", las=1, font.lab=2)
