library(cyclar)

# Sliding window functions

# Generate time series

x <- arima.sim(n=100, model=list(ar=c(0,-0.7)))
x <- ts(x, start=1910)

x <- ar2.sim(n=100, mu=2, a=-1, a2=-0.7)
ts.summary(x, plot=TRUE)

x.n <- ts(replicate(arima.sim(n=100, model=list(ar=c(0,-0.7))), n=4))

# For single time series
ts.window(x=x, width=16, fun="mean", plot.type="single")
ts.window(x=x, width=16, fun="mean", plot.type="multiple")

# For multiple time series
ts.window.n(x=x.n, width=16, fun="mean", plot.type="single")

# Multiple series, plot.type="multiple"
x.m <- ts.window.n(x=x.n, width=16, fun="mean", plot.type="multiple")
x.s <- ts.window.n(x=x.n, width=16, fun="sd", plot.type="multiple")

# ACF
p <- ar2.parms(k=5,v=6,output=list)
n <- 200
x <- arima.sim(n=n, model=list(ar=c(1+p$a1,p$a2)), sd=.3) + 2
x <- ts(x, start=1900)

par(mfrow=c(2,1))
	plot(x)
	plot(exp(x))
par(mfrow=c(1,1))

ar2.plot(k=5,v=6)
points(1+p$a1, p$a2, pch=16, col=2)

arima(x,order=c(2,0,0))
ar.yw(x)

# Create custom function that returns the acf-function
my.acf <- function(x,lag.max=NULL) {
	out <- acf(x, plot=FALSE, lag.max=lag.max)
	z <- out$acf[,,1]
	names(z) <- out$lag[,,1]
	z
}

my.acf(x)
# Use rollapply
y <- rollapply(as.zoo(x), width=17, FUN="my.acf", lag.max=16)
t(y)

# Use ts-window, the wrapper-function for rollapply
f <- ts.window(x, width=17, fun="my.acf", lag.max=16, plot=FALSE) # plot doesn't work with this function
t(f$window)

ts.diag.acf(x)

# Extract acf for period length 4 & 5
my.acf.45 <- function(x,lag.max=NULL) {
	z <- my.acf(x,lag.max)[c("4","5")]
	z
}

my.acf.45(x)

z1 <- rollapply(as.zoo(x), width=17, FUN="my.acf.45", lag.max=16)
matplot(z1, type="l", lwd=2)

z <- ts.window(x, width=16, fun="my.acf.45", plot=FALSE)

plot(x = range(time(z$window)), y = c(0,1), type="n", xlab="Time", ylab="ACF", font.lab=2, las=1)
points(time(z$window), z$window[,1], type="b", col=1, pch=16)
points(time(z$window), z$window[,2], type="b", col=2, pch=1)
