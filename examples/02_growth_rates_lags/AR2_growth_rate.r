################################################################################
library(cyclar)
################################################################################

# Generate data from an AR(2)-process

# Select a series with negative direct density dependence, i.e. (1+a1) < 1, a1 < 0
k <- 3
v <- 4

p <- ar2.parms(k=k, v=v)
ar2.period(p$a1,p$a2)
ar2.plot(k=k,v=v)
points(1+p$a1, p$a2, pch=16, col=2)

# Simulate one time series
a0 <- 0; a1 <- p$a1; a2 <- p$a2
ar.var <- 0.2
n <- 100

# x.sim <- ar2.sim(n=n, a0=a0, a1=a1, a2=a2, sd=sqrt(ar.var)); x <- x.sim$ts.x
x <- arima.sim(n=n, model=list(ar=c(1+a1,a2)), sd=sqrt(ar.var))
x <- ts(x,start=1910)

range(x)
ts.summary(x,plot=TRUE)
ts.spec(x)

# Growth rate functions
# ar2.lag, ar2.lag.plot, ar2.gr

# Create data frames with lagged variables
# ar2.lag assumes that data is on log-scale!

N <- exp(x) # natural scale
z1 <- ar2.lag(x=x, log=FALSE)
z2 <- ar2.lag(x=N, log=TRUE) # log-transform the series N
all.equal(z1,z2,tol=0)

# Plot growth rates on "natural" scale
ar2.gr(x=N)
ar2.gr(x=lynx, method="lambda")

par(mfrow=c(1,2))
hist(lynx, col="steelblue", breaks=20,freq=FALSE,xlab="Lynx",main="Mackenzie lynx series")
lines(density(lynx),col=2)
hist(log(lynx),col="steelblue", breaks=20,freq=FALSE,xlab="log(Lynx)",main="Mackenzie lynx series")
lines(density(log(lynx)), col=2)
par(mfrow=c(1,1))

# Lag plots
ar2.lag.plot(x, lag.max=5, variable="r")
ar2.lag.plot(x, lag.max=5, variable="x")

# Example with the lynx data set
ar2.lag.plot(x=lynx, variable="x", lag.max=8, log=TRUE, base=10, constant=0, pch=16, cex=0.8)
ar2.lag.plot(x=lynx, variable="r", lag.max=8, log=TRUE, base=10, constant=0, pch=16, cex=0.8)

lag.plot(as.numeric(log10(lynx)), lags=8) # compare with

ar2.lag.plot(x=lynx, variable="r", lag.max=3, log=TRUE, base=10, constant=0, pch=16, cex=0.8)
ar2.lag.plot(x=lynx, variable="x", lag.max=4, log=TRUE, base=10, constant=0, pch=16, cex=0.8)

# Compare
# lag.plot uses positive lags, which in fact are "leads" rather than "lags"
z <- ar2.lag.plot1(x=lynx, variable="x", lag=1, log=TRUE, base=10, constant=0, pch=16, cex=0.8)
lag.plot(log10(lynx),lag=1,do.lines=FALSE)
points(z$y,z$x,col=2,pch="+") # x & y switch because lag.plot uses postive lags, "leads"

################################################################################

# 3-dimensional of 'reproduction surfaces'
r.3d <- function(x,y,a0,a1,a2) a0 + a1*x + a2*y

# Persp plot
res <- curve3d(r.3d(x=x,y=y,a0=1,a1=a1,a2=a2),phi=30,theta=60,xlab=expression(t[-1]),ylab=expression(t[-2]),
	zlab="R",ticktype="detailed",col="steelblue")

# contours (not very interesting for linear surfaces)
#curve3d(r.3d(x=x,y=y,a0=1,a1=a1,a2=a2),xlab=expression(t[-1]),ylab=expression(t[-2]),zlab="R",sys3d="contour")

################################################################################

# TEST SECTION - not complete!

s.lynx <- ts.stand(log(lynx))
plot(s.lynx)
ts.summary(s.lynx)
ar2.lag.plot(s.lynx,lag.max=8)

# Can ar2.lag be used with scaled, standardized data?
z <- ar2.lag(s.lynx,lag.max=4)

fm1 <- arima(s.lynx,order=c(2,0,0))

x <- seq(-3, 3, length=50)
y <- x
f <- function(x,y,a0,a1,a2) a0 + a1*x + a2*y
z <- outer(x, y, f, a0=0, a1=coef(fm1)[1]-1, a2=coef(fm1)[2])

op <- par(bg = "white")
res <- persp(x, y, z, theta = -50, phi = 30, col = "lightblue", ticktype="detailed",
	 xlab = "X", ylab = "Y", zlab = "Growth rate")

round(res, 3)

# Add points
points(trans3d(z[,2], z[,3], z[,6], pmat = res), col = 2, pch =16)

################################################################################
