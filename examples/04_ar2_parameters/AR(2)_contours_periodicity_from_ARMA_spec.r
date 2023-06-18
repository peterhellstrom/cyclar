################################################################################
library(cyclar)
################################################################################

# Load Shumway & Stoffer's workspace:
load("C:\\WORK\\ANALYSER\\-= R Development =-\\cyclar\\Examples\\tsa3.rda")

# The expected periodicity from the ARMA spectrum is NOT the same as the period contours given by Royama.

# I have used:
# Shumway & Stoffer's workspace tsa3.rda and the function spec.arma.
a1 <- 0
a2 <- -0.7
z <- spec.arma(ar=c(a1,a2), n.freq=1001)
z$freq[which(z$spec == max(z$spec))]

# I changed the original function a little bit (removed ma, & error trapping. I also created a slightly different plot.)
# The new function is called ts.spec.ar2.

# Example 1
# Create parameters for a 5-year cycle
k <- 5
v <- 1.25
p <- ar2.parms(k=k, v=v, plot=FALSE, output=list, kmethod="royama")

# But notice that peak of spectral density is not located at 1/k (0.2)
# The period obtained from the spectral density is slightly longer than 5.
x1 <- ts.spec.ar2(a1=p$a1,a2=p$a2, sd=1, plot=TRUE)
x1$f; x1$period

# Plot both Royama & Jenkins/Watts contours to illustrate this.
ar2.plot(k=NULL,v=v)
ar2.k.add(k, col=2)
ar2.k.add(k=x1$period, kmethod="jenkins", col=3, lty=2)
ar2.k.add(k=k, kmethod="jenkins", col=4, lty=2)
with(p, points(1+a1,a2,pch=16,col=2))

# Check simulation
x2 <- arima.sim(n=1000, model=list(ar=c(p$phi1,p$phi2)))
spec.ar(x2, order=2) # Check estimated spectral density
# Add theoretical expectation in red
points(x1$freq, x1$spec, type="l",col=2,lty=3)
abline(v=x1$f, col=2, lty=2)
legend("topright", c("Simulated","Theoretical"), col=c(1,2), lty=c(1,3), bty="n")


# Example 2
# If variance increases, then Royama's contour is almost identical to peak of spectral density.
# 5-year cycle with ipv=4
a1 <- 0.5700546
a2 <- -0.8507622
ar2.period(a1,a2)
1/ar2.period(a1,a2)[,"period"]

x1 <- ts.spec.ar2(a1=a1,a2=a2, plot=TRUE)
x1$f; x1$period

# Note that k-contours are very similar at high variances
ar2.plot(k=NULL, v=ar2.ipv(a1,a2))
ar2.k.add(k=ar2.period(a1,a2)[,"period"],lty=1,col=2)
ar2.k.add(k=x1$period,lty=2,col=3,kmethod="jenkins")
points(1+a1,a2,pch=16,col=2)

# Find contours for quasi-period in AR(2)-space:
# First a numerical evaluation, then an exact solution
library(emdbook)
lev <- seq(0,10,0.5)

# This of course evaluates the function in a lot of impossible regions, but ignore warnings...
f1 <- curve3d(ts.spec.ar2(a1=x-1,a2=y,plot=FALSE)$period, from=c(-2,-1), to=c(2,0), n=c(100,100), sys3d="contour", levels=lev)

# Plot in AR(2)-space
k <- c(2,3,3.5,4,4.5,5,5.5,6,8)
ar2.plot(k=k, v=5, vlab=TRUE, vcontours=TRUE, vcol="darkgrey", cex.contours=0.8)
contour(f1, levels=k, col=4, lty=2, add=TRUE, cex=1)

# The analytical solution:
# Cryer & Chan (2008), p. 336 - 338, present the work of Jenkins and Watts (1968)
# See especially equation 13.5.9, p. 338

# Take equation 13.5.9 and solve for a2 and we up with:
f2 <- function(a1,k) -(1+a1) / (4*cos(2*pi/k) - (1+a1))

# This has been incorporated in ar2.plot
ar2.plot(kmethod="royama")
ar2.plot(kmethod="jenkins")

# add contours to existing plot
kj <- c(3,5,6,8)
ar2.plot(k=NULL,v=NULL)
for (i in 1:length(kj)) ar2.k.add(k=kj[i], col=3, kmethod="jenkins", lty=2)
for (i in 1:length(kj)) ar2.k.add(k=kj[i], col=2, kmethod="royama")

curve(f2(a1=x-1,k=5),add=TRUE,col=5)
