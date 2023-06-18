# Calculate moving averages (and other statistics) with ts.moving:
library(cyclar)

x <- ts(runif(n=16),start=1980)

ts.moving(x, output="vector", names="span")
ts.moving(x, output="ts")

ts.moving(x,width=5)
# compare with zoo-function rollmean()
rollmean(x, k=5)

# Test another function, standard deviation
ts.moving(x,width=5,fun="sd")

# Plot moving averages of different window widths (3,5,7 & 15 points):
x <- ts(runif(n=100,0,10))

plot(x, pch=16, type="p", main="Moving average")
lines(x, lty=2)
lines(ts.moving(x,width=3,output="ts"),col=2,lwd=2)
lines(ts.moving(x,width=5,output="ts"),col=3,lwd=2)
lines(ts.moving(x,width=7,output="ts"),col=4,lwd=2)
lines(ts.moving(x,width=15,output="ts"),col=5,lwd=2)
legend("topleft", legend=c(3,5,7,15), lty=c(1,1,1,1), lwd=c(2,2,2,2), col=c(2,3,4,5), bty="o", title="window", cex=0.7, bg="white")

# Plot ACF of series x, and three-point moving average
# Note that for a series of randomly drawn (uniform) numbers, moving averages with few data points can introduce spurious
# periodicity in the series of moving averges. Check acf's:
par(mfrow=c(1,2))
acf(x)
acf(ts.moving(x,width=3,output="ts"))
par(mfrow=c(1,1))

# Simulate a nonstationary AR(1)-series
x <- arima.sim(n=100, model=list(ar=0.9)) # Random walk (non-stationary)
x <- ts(x, start=1910)
x.ma <- ts.moving(x, na.rm=TRUE, width=10, fun="mean", output="ts")
x.ma

plot(x,lty=2, main="Moving average of a random walk")
points(x, pch=16, col=1, cex=0.7)
lines(x.ma, lty=1, col=2, lwd=2)

# Simulate AR(2) & ARMA(2,1)
x.ar <- arima.sim(n=200, model=list(ar=c(0,-0.7)))
x.arma <- arima.sim(n=200, model=list(ar=c(0,-0.7), ma=c(0.9)))
x.ar.sd <- ts.moving(x.ar, na.rm=TRUE, width=13, fun="sd", output="ts")
x.arma.sd <- ts.moving(x.arma, na.rm=TRUE, width=13, fun="sd", output="ts")

par(mfrow=c(2,1))
plot(x.ar, lty=2, col=4, main="Moving SD of an AR(2)")
lines(x.ar.sd, lty=1, col=4, lwd=2)
plot(x.arma, lty=2, col=2, main="Moving SD of an ARMA(2,1)")
lines(x.arma.sd, lty=1, col=2, lwd=2)
par(mfrow=c(1,1))

# Handle missing values (NAs)
x[c(5,15)] <- NA
x

ts.moving(x)
ts.moving(x, na.rm=TRUE)

# Handle vectors (i.e. not time-series)

x <- runif(n=100)
ts.moving(x)
