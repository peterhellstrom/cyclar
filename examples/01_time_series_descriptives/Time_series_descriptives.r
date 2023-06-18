library(cyclar)

# Test the (new versions of the) time series descriptive functions

# Generate a random time series, just uniform numbers
#N <- ts(runif(100000,min=3,max=5))
N <- ts(runif(100,min=3,max=5))

N <- arima.sim(n=100, model=list(ar=c(0,-0.7)))

par(mfrow=c(1,2))
	hist(N, breaks=30, col="steelblue")
	hist(diff(N), breaks=30, col="steelblue")
par(mfrow=c(1,1))

ts.summary(N)
ts.stand(N)

ts.detr(N)
# ts.detr.step(N, g, k=)

ts.cobweb(N)
ts.binary(N)

ts.spec(N)
ts.spec(diff(N))

ts.diag.acf(N, f="acf", max.only=FALSE, output=rbind)
ts.diag.acf(N, f="pacf", max.only=TRUE, plot=TRUE)
ts.diag.acf(N, fun="pacf", max.only=TRUE)

par(mfrow=c(1,2))
	ts.diag.acf(N,f="acf",max.only=TRUE,plot=TRUE)
	ts.diag.acf(N,f="pacf",max.only=TRUE,plot=TRUE)
par(mfrow=c(1,1))

par(mfrow=c(1,2))
	ts.diag.acf(N,plot=TRUE)
	ts.diag.acf(diff(N),plot=TRUE)
par(mfrow=c(1,1))

par(mfrow=c(1,2))
	ts.diag.acf(N, f="pacf", plot=TRUE)
	ts.diag.acf(diff(N), f="pacf", plot=TRUE)
par(mfrow=c(1,1))
