library(cyclar)

# Sinus curve
xv <- seq(0, 100, 0.1)
yv <- sin(xv)
yv <- ts(yv, start=1978.0, deltat=1/12)

yv <- ts(yv, start=1978.0, deltat=1)

plot.ts(yv, bty="l", xlab="Time", ylab=expression(sin(x)), font.lab=2, main="Sinus curve")

length(yv)
acf.yv <- acf(yv, lag.max=500, plot=TRUE)
str(acf.yv)

# Use the function ts.diag.acf
ts.diag.acf(yv, lag.max=40)
ts.diag.acf(yv, lag.max=120)
ts.diag.acf(yv, lag.max=500)

ts.diag.acf(yv, lag.max=500, output=rbind)
z <- ts.diag.acf(yv, lag.max=500, output=cbind)

plot(acf.yv$lag, acf.yv$acf, type="l", xlab="Lag", ylab="Autocorrelation coefficient", font.lab=2, bty="l", main="Find dominant period")
abline(h=0, lty=2, col=1)
points(z[1,"max"],z[2,"max"],col=2,pch=16)
points(z[1,"min"],z[2,"min"],col=4,pch=16)

# Generate another example
x <- arima.sim(n=100, model=list(ar=c(0,-0.7)))
# x <- ts(runif(n=100,10,20))

ts.spec(x)
ts.diag.acf(x)
ts.diag.acf(x,fun="pacf")
