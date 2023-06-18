################################################################################
library(cyclar)
################################################################################

# Check/draw ACF for a stationary AR(2)-process in the region with complex roots
k <- 5; v <- 4
p <- ar2.parms(k=k, v=v, output=list)
ar2.acf(a1=p$a1, a2=p$a2)

ar2.acf(lags=0:10, a1=p$a1, a2=p$a2, plot=TRUE)
ar2.acf(a1=p$a1, a2=p$a2, plot=TRUE, subplot=FALSE)
ar2.acf(a1=p$a1, a2=p$a2, plot=TRUE, subplot=TRUE)

ar2.acf(a1=p$a1, a2=p$a2, plot=TRUE, subplot=TRUE, hadj=-2)

# Compare with ARMAacf.
lag.max <- 40
ar2.acf(lags=0:lag.max, a1=p$a1, a2=p$a2, plot=FALSE)

ar2.acf(a1=p$a1, a2=p$a2, plot=TRUE,lag.max=lag.max)
ARMAacf(ar=c(1+p$a1,p$a2),lag=lag.max)

plot(0:lag.max, ARMAacf(ar=c(1+p$a1,p$a2),lag=lag.max), type="b")
points(0:lag.max, ar2.acf(lags=0:lag.max, a1=p$a1, a2=p$a2, plot=FALSE), type="b", col=2, lty=2)

library(rootSolve)

f.acf <- function(x,phi1,phi2) ((sqrt(-phi2)^x)*(sin(acos(phi1/(2*sqrt(-phi2)))*x + atan(((1-phi2)/(1+phi2)))) / sin(atan(((1-phi2)/(1+phi2))))))

curve(f.acf(x=x,phi1=0,phi2=-0.7), from=0, to=40, n=1001, xlab="lag", ylab="acf")

parms <- list(phi1=0, phi2=-0.7)
z <- fderiv(f=f.acf, parms=parms, tangent=FALSE)



phi1 <- 0
phi2 <- -0.7
x <- arima.sim(n=100, model=list(ar=c(phi1,phi2)))

# http://en.wikipedia.org/wiki/Autocorrelation
# Fast Fourier Transform
fourier <- fft(x)
# extract the power spectrum (power is sometimes referred to as "magnitude")
# magnitude=sqrt(Re(fourier)*Re(fourier)+Im(fourier)*Im(fourier))
magnitude <- Mod(fourier)
# extract the phase which is atan(Im(fourier)/Re(fourier))
phase <- Arg(fourier)
