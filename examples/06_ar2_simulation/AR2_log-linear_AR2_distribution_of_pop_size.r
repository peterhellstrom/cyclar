library(cyclar)

# Q: How much variability in mean & median should we expect
# when simulating an AR(2)-process a large number of times:

# Start with setting up one simulation
n <- 10000
mu <- 2.5 # mean on log-scale
# Generate parameters
p <- ar2.parms(k=3, v=2, output=list)
a1 <- p$a1
a2 <- p$a2
# Calculate intrinsic process variance
ipv <- ar2.ipv(a1,a2)
# Set white noise of innovations
v <- 0.3

# Generate one sample
x <- arima.sim(n=n, model=list(ar=c(1+a1,a2)), sd=sqrt(v)) + mu
range(x)

# Compare empirical and theoretical mean:
mean(x)
mu

# Compare empirical and theoretical sd:
sd(x)
sqrt(ipv * v)

# Plot distribution of log population size (should appear normal):
hist(x, col="steelblue", breaks=30, freq=FALSE)
lines(density(x)) # add kernel
curve(dnorm(x, mean=mu, sd=sqrt(ipv * v)), col=2, lwd=2, add=TRUE) # theoretical distribution
abline(v=mu, col=2, lwd=2, lty=2) # expected mean
abline(v=mean(x)) # observed mean


# Exponentiated series

# Compare empirical and theoretical mean:
mean(exp(x))
exp(mu + 0.5*(v*ipv))

# Compare empirical and theoretical sd:
sd(exp(x))
sqrt(exp(2*mu + (v*ipv))*(exp(v*ipv)-1))

# Repeat the simulations a large number of times

nrepl <- 1000 # Replicate simulation this number of times
xn <- replicate(arima.sim(n=n, model=list(ar=c(1+a1,a2)), sd=sqrt(v)) + mu, n=nrepl)

# Plot distribution of difference from expected mean
xn.diff <- colMeans(xn) - mu
range(xn.diff)
hist(xn.diff, col="steelblue", breaks=30)

# Exponentiated series
xn.exp <- exp(xn)

# mean
xn.exp.diff <- colMeans(xn.exp) - exp(mu + 0.5*(v*ipv))
range(xn.exp.diff)
hist(xn.exp.diff, col="steelblue", breaks=30)

# median
# For exp(x), note that the mean is exp(mu)
hist(apply(xn.exp,2,median), col="steelblue", breaks=30)
abline(v=exp(mu), col=2, lty=2, lwd=2)

xn.exp.median.diff <- apply(xn.exp,2,median) - exp(mu)
range(xn.exp.median.diff)
hist(xn.exp.median.diff, col="steelblue", breaks=30)
