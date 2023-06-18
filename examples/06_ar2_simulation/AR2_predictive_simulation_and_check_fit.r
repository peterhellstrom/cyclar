################################################################################
library(cyclar)
################################################################################
# Code borrowed from Gelman & Hill, p. 163-165, then updated:
################################################################################

# Using predictive simulation to check the fit of a time-series model

# Generate a time series
k <- 5
v <- 3

p <- ar2.parms(k=k,v=v)
ar2.period(a1=p$a1,a2=p$a2)

mu <- 2
a1 <- p$a1
a2 <- p$a2
ar.var <- 0.2
n <- 100

x <- arima.sim(list(order=c(2,0,0), ar = c(1+a1,a2)), n=n, sd=sqrt(ar.var)) + mu

# Test some other series, pick one series and re-run code
x <- lynx
x <- kilpis[,"cr.f"]

# Estimate the coefficients with the Yule-Walker equations.
fm1 <- ar.yw(x, order.max=2)

################################################################################

# Use estimated coefficients as inputs for simulations
b.hat <- fm1$ar # vector with coefficients
s.hat <- sqrt(fm1$var.pred) # standard deviation
x.mean <- fm1$x.mean # mean
n <- length(x) # length of each series
n.sims <- 10000 # number of replicate simulations

# Simulate replicated datasets
# Create a function
sim <- function(x,n,mu,a1,a2,sd) {

	y <- numeric(n)
	p <- numeric(n)
	y[1:2] <- x[1:2]

	b.hat <- c(1+a1,a2)
	s.hat <- sd

	for (t in 3:n) {
		prediction <- (c(y[t-1], y[t-2]) %*% b.hat) + mu
		y[t] <- rnorm(1, prediction, s.hat)
	}
	as.ts(y)
}
################################################################################
# Do the simulations
system.time(z <- replicate(sim(x=x,n=n,mu=x.mean,a1=b.hat[1]-1,a2=b.hat[2],sd=s.hat),n=n.sims))

# Check a subset of randomly selected time series
inds <- sample(1:n.sims, 25, replace=F)
z.plot <- z[,inds]

par(mfrow=c(5,5))
op <- par(mar=c(1,2,1,1))
for (i in 1:25) plot(z.plot[,i],type="l",xlab="",ylab="",xaxt="n")
par(mfrow=c(1,1))
par(op)

# Construct test statistic that is the frequency of switches.
# The number of years in which an increase is followed by a decrease, or vice versa.

Test <- function(y) {
	n <- length(y)
	y.lag <- c(NA, y[1:(n-1)])
	y.lag2 <- c(NA, NA, y[1:(n-2)])
	sum(sign(y - y.lag) != sign(y - y.lag2), na.rm=TRUE)
}

# Calculate this test-statistic for each series in matrix z
tr <- apply(z,2,Test)
# Plot distribution of test-statistic
plot(density(tr), main="Distribution of test statistic")
# Add line for input series x
abline(v=Test(x), col=2, lty=2, lwd=2)

# The input series x has this number of switches:
Test(x)

# Evaluate how many of the simulated series that have more switches than series x
table(tr > Test(x)) / n.sims

# For the lynx series, ~84% of the simulated series had more switches than the lynx series.
# For the Kilpisj?rvi data series for fall densities of grey-sided voles, ~95% of the
# simulated series has fewer switches than the observed data, suggesting that the AR(2) model has a poor fit.
