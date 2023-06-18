################################################################################
library(cyclar)
################################################################################

p <- ar2.parms(k=3, v=3)
mu <- 2.1
a1 <- p$a1 # direct d-d
a2 <- p$a2 # delayed d-d
ar2.var <- 0.3 # white noise

n <- 10000

# Demonstrate the two functions available in cyclar, ar2.sim & ar2.sim2 (which can take innovations in vector form as an argument) and arima.sim:

# with ar2.sim
z1 <- ar2.sim(mu=mu,a1=a1,a2=a2,n=n,sd=sqrt(ar2.var))
ts.summary(z1)
sqrt(ar2.ipv(a1,a2,sd=sqrt(ar2.var)))

# with ar2.sim2
z2 <- ar2.sim2(mu=mu,a1=a1,a2=a2,n=n,sd=sqrt(ar2.var))
ts.summary(z2)
sqrt(ar2.ipv(a1,a2,sd=sqrt(ar2.var)))

################################################################################
# various calls to arima.sim

# Simulation of an AR(2) process is usually done with the arima.sim function,
# which calls the filter function, and then filters the innovations with the ar-weights in a recursive fashion.
# Different calls to arima.sim for generating an AR(2)-process:

# if not supplied, a reasonable value for n.start is calculated internally.
x <- arima.sim(list(order=c(2,0,0), ar = c(a1,a2)), n=n, sd=sqrt(ar2.var), n.start=500, rand.gen=rnorm) + mu
# With innovations generated before call to arima.sim:
w <- rnorm(n=n, sd=sqrt(ar2.var))
x1 <- arima.sim(list(order=c(2,0,0), ar = c(a1,a2)), n=n, innov=w, n.start=500) + mu
# The component order is optional:
x2 <- arima.sim(model=list(ar = c(a1,a2)), n=n, sd=sqrt(ar2.var), n.start=500) + mu

################################################################################
# in ar2.sim2, method="filter" is equal to method="recursive".
# To show this, we can set up two simulations, one with each method, and supply the same
# start innovations and innovations. It is simportant to supply BOTH start.innov and innov.
# Supplying only innov would generate two different random sequences of start.innov, and
# the time series would differ.

mu <- 0
a1 <- -0.5
a2 <- -0.7

n <- 1000
n.start <- 500
w <- rnorm(n=n)
start.w <- rnorm(n=n.start)

zf <- ar2.sim2(mu=mu, a1=a1, a2=a2, n=n, innov=w, start.innov=start.w, method="filter")
zr <- ar2.sim2(mu=mu, a1=a1, a2=a2, n=n, innov=w, start.innov=start.w, method="recursive")

# Summary for zf and zr should be identical
ts.summary(zf)
ts.summary(zr)

# Plot difference between zf and zr
plot(zoo(zf-zr))
all.equal(zf,zr,tol=0)

head(zf + 2,50)
head(zr + 2,50)

# ... to be cont.
################################################################################
# ... cont. from previous section
# This section shows that the filter and recursive method are equivivalent
# when generating a log-linear AR(2)-process.

mu <- 2.1
a1 <- -1.25
a2 <- -0.7
n <- 100

# Generate innovations
z <- rnorm(n=n)

# Filter innovations
y <- filter(z, c(1+a1,a2), method="recursive") + mu

# Compare with recursive equation
# This calculation should give exactly the same result as the filter-method.
# Note that different expressions are necessary for t[1] & t[2].
# I also tried ignoring the t==2 part in the loop. If t==2 is left out,
# the two series converge after roughly n.start steps (see calculation of n.start below).

# similar code is used in ar2.sim2:
a0 <- mu*(-(a1+a2))
x <- z + mu #  generate vector x, and the first value of this vector is actuall y kept, all other are recursively replaced.
for (t in 2:length(z)) {
	if (t == 2) x[t] <- mu*(-a1) + (1+a1)*(x[t-1]) + z[t] # Note that the value for t=2 is actually AR(1)
	else x[t] <- a0 + (1+a1)*x[t-1] + a2*x[t-2] + z[t]
}

# Calculate n.start
ar <- c(1+a1,a2)
p <- length(ar)
minroots <- min(Mod(polyroot(c(1, -ar))))
n.start <- p + ifelse(p > 0, ceiling(6/log(minroots)), 0)

# Plot time series y [filter] & x [recursive]
plot(y,type="n",xlab="Time",ylab="Simulated series, x & y")
lines(y)
lines(x,col=2,lty=2)
legend("topleft",c("y [filter]","x [recursive]"), col=1:2, lty=1:2, bty="n", cex=0.8)
abline(v=n.start, col=4, lty=3)

# Plot difference between y & x
plot(zoo(y-x),ylab="y-x")
abline(v=n.start, col=4, lty=3)
all.equal(y,as.ts(x),tol=0)

# Relationship between a0 (intercept) and mu (mean).
# Show an example, plot a0 against (1+a1) and a2:

curve3d(ar2.intercept(mu=2, a1=x-1, a2=y), from=c(-2,-1), to=c(2,1), n=c(20,10), col="steelblue",
ticktype="detailed", theta=50,zlab="a0",xlab="1+a1",ylab="a2")
title(main=bquote(paste("AR(2) intercept [", a[0], "] for ", mu == .(mu),sep="")))

ar2.intercept(mu=2, a1=-2, a2=-0.5)
ar2.intercept(a0=5, a1=-2, a2=-0.5)

x <- ar2.sim2(n=100, mu=2, a1=-2, a2=-0.5)
attr(x,"p")

# Innovations in arima.sim
# Innovation are generated in the following in the following way in arima.sim
# and arima.sim2.
# In the example below, the innovations are drawn from a normal distribution with mean = 0 & sd = 1.
# Two sets of innovations are created, first innovations for the burn-in sequence (start.innov)
# and then innovations starting at t=1. The burn-in sequence is discarded before output is printed.

n.start <- 50
n <- 100
start.innov <- rnorm(n.start)
innov <- rnorm(n)
# Bind start.innov & innovations together (the arima.sim function then filters this ts).
x <- ts(c(start.innov, innov), start = 1 - n.start)

# After filtering, the "burn-in"-sequence is removed:
as.ts(x[-(1L:n.start)])
# If x is a ts-object, this could also be achieved with:
window(x,1,100)
