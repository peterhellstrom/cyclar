##################################################################
library(cyclar)
library(TSA)
##################################################################

# Periodograms

# Generate a time series, AR(2)
k <- 6
v <- 3
p <- ar2.parms(k=k,v=v,plot=TRUE)
a1 <- p$a1
a2 <- p$a2
mu <- 2

n <- 100
x <- ar2.sim(n=n, mu=mu, a1=a1, a2=a2, sd=sqrt(0.2))
ts.summary(x,plot=TRUE)
ts.spec(x)

# Periodogram

# Functions spec.pgram & spectrum
# spectrum is a wrapper that calls spec.pgram

?spec.pgram

# Vary the degree of smoothing of the periodogram by changing the parameter "spans". Provide odd integers?
# Try to vary detrend=T and demean=F 
# spans = width of Daniell smoothers

# most common arguments arguments
# pad = proportion of data to pad. Zeros are added to the end of the series to increase its length by the proportion pad.
# Defaualt is 0
# spans = vector of odd integers giving the width of Daniell smothers. Default is NULL
# log = 'no'
tmp1 <- spec.pgram(x, spans=c(3,3))
str(tmp1)

# Test the wrapper function spectrum
?spectrum

tmp2 <- spectrum(x) # additional arguments can be passed, see spec.pgram
str(tmp2)

# How to find the peak frequencies?
# The function "peaks" can be used to find a local maxima.
# For objects of class 'spec', 

# Find peak frequencies
# default span width in ts.spec.peaks is 9
ts.spec.peaks(tmp1, span=3)
ts.spec.peaks(tmp1)

ts.spec.peaks(tmp2)

# Bivariate example, spec.pgram
xmat <- ar2.sim.npatch(n=n, npatch=2, mu=mu, a1=a1, a2=a2, rho=0.9, sd=1)$x
tmpbv <- spec.pgram(ts.union(xmat[,1],xmat[,2]), spans=c(3,3))

plot(tmpbv, plot.type="coherency")
plot(tmpbv, plot.type="phase")




# Create a function that generates a series x, then estimates the spectrum and return the peak periodicity:
fn <- function() {
	x <- ar2.sim(n=n, mu=mu, a1=a1, a2=a2)
	tmp <- spec.pgram(x, spans=c(3,3),plot=FALSE)
	as.numeric(ts.spec.peaks(tmp, span=5, plot=FALSE)[1,2])
}

# Test the function once
fn()

# Replicate a 'large' number of times and view distribution of estimated periods
z <- replicate(fn(), n=1000)
summary(z)
hist(z, breaks=30, col="steelblue", freq=FALSE)
# Expected line
abline(v=1/(acos(as.numeric(-((1+a1)*(1-a2))/(4*a2)))/(2*pi)), col=2, lwd=2, lty=2)


# in TSA
periodogram(x)

# Fast Fourier transform
per <- abs(fft(x))^2


# Compare theoretical with empirical spectral density

# Generate a new series x, where innovations are z ~ N(0,1)
# (otherwise, the empirical and theoretical densities do not match).
x <- ar2.sim(n=n, mu=0, a1=a1, a2=a2)

phi <- c(1+a1,a2)
sp <- spec(x,log='no',xlab='Frequency', ylab='Sample Spectral Density',sub='',plot=F)
plot(sp)
lines(sp$freq, ARMAspec(model=list(ar=phi),freq=sp$freq,plot=F)$spec, lty='dotted')
abline(h=0)


k <- kernel('daniell', m=5)
sp <- spec(x,kernel=k,log='no',sub='',xlab='Frequency', ylab='Smoothed Sample Spectral Density', plot=F)
plot(sp)
lines(sp$freq,ARMAspec(model=list(ar=phi),freq=sp$freq, plot=F)$spec,lty='dotted')


sp1 <- spec(x,spans=3,sub='',lty='dotted',xlab='Frequency', ylab='Log(Estimated Spectral Density)')
sp2 <- spec(x,spans=9,plot=F)
sp3 <- spec(x,spans=15,plot=F)
lines(sp2$freq,sp2$spec,lty='dashed',col=2)
lines(sp3$freq,sp3$spec,lty='dotdash',col=4)
f <- seq(0.001,.5,by=.001)
lines(f,ARMAspec(model=list(ar=c(1+a1,a2)),freq=f, plot=F)$spec,lty='solid')
