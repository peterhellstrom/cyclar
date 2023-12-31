library(cyclar)
library(pastecs)
library(TeachingDemos)
library(zoo)
# Q: Is it possible that shorter cycles 3-4 years are (falsely) detected as 7-year cycles?

# 1) Generate some vectors with states 0,1 or 0,1,2,3
# 2) check corresponding acf's and estimate AR-paramaters from acf (acf2AR).
# 3) Calculate AR(2)-dynamics for parameters estimated in step 2.
# 4) Optional, plot theoretical acf for the AR(2)-process generated by parameters in step 2/3.

# Length of series
n <- 30

# Sequence a
x1 <- rep(c(c(0,0,1),c(0,0,1,1)),length.out=n)
p1 <- acf2AR(acf(x1, plot=FALSE)$acf)
phi1 <- p1[2,1:2]
ar2.period(a1=phi1[1]-1, a2=phi1[2])
z1 <- ar2.acf(a1=phi1[1]-1, a2=phi1[2], plot=TRUE)
ar.acf(ar=phi1)
points(as.numeric(names(z1)), z1, col=2)

# Sequence b
x2 <- rep(c(c(0,1,3),c(0,1,2,3)),length.out=n)
p2 <- acf2AR(acf(x2, plot=FALSE)$acf)
phi2 <- p2[2,1:2]
ar2.period(a1=phi2[1]-1, a2=phi2[2])
ts.binary(x2)==x1

# Sequence c
x3 <- rep(c(0,1,2,3,2,1), length.out=n)
p3 <- acf2AR(acf(x3, plot=FALSE)$acf)
phi3 <- p3[2,1:2]

# Plot all series
par(mfrow=c(3,2))
	plot(x1,type="b",pch=16,main="Sequence 3a",xlab="Time",ylab="Index",font.lab=2,las=1,bty="l")
	acf(x1)
	plot(x2,type="b",pch=16,main="Sequence 3b",xlab="Time",ylab="Index",font.lab=2,las=1,bty="l")
	acf(x2)
	plot(x3,type="b",pch=16,main="Sequence 3b",xlab="Time",ylab="Index",font.lab=2,las=1,bty="l")
	acf(x3)
par(mfrow=c(1,1))

# Lag plots reveal why lag=7 is dominant for series a & b, and 3 & 6 for series c.
lag.plot(x1, lags=12)
lag.plot(x2, lags=12)
lag.plot(x3, lags=12)

# Markov-type estimation of cycle length gives ~correct results
state.cycle(phase(x1))
state.cycle(phase(x2))
state.cycle(phase(x3))

ts.summary(x1, plot=TRUE)
ts.summary(x2, plot=TRUE)
ts.summary(x3, plot=TRUE)

ts.spec(ts(x1))
ts.spec(ts(x2)) # note bimodal spectral density
ts.spec(ts(x3))

# Plot AR(2) coefficients
ar2.plot(ylim=c(-1,0), triangle=FALSE, width=9, height=5, text.lab=FALSE)
points(phi1[1], phi1[2], col=2, pch=16, cex=1.3)
points(phi2[1], phi2[2], col=2, pch=16, cex=1.3)
text(x=phi1[1], y=phi1[2]-0.05, labels="a", font=2, cex=1.2)
text(x=phi2[1], y=phi2[2]-0.05, labels="b", font=2, cex=1.2)


################################################################################
################################################################################
# ACF BIFURCATIONS

# The problem:
# If cycle length is non-integer, such as 2.5 or 3.5, using the autocorrelation function to find cycle length
# can give the unexpected result that the maximum (dominant) peak of the acf is located at double the true period.
# This occurs in certain regions (following a bifurcating pattern) in the a1/a2 space of an AR(2)-process.
# These regions are characterized by high intrinsic process variance, where the acf does not dampen quickly.

# A numerical example, 3.5 year cycle with very high ipv:
k <- 3.5
v  <- 20
p <- ar2.parms(k=k, v=v, plot=TRUE)

# Check theoretical acf, notice that it appears that we have a 7-year cycle.
z1 <- ar.acf(ar=c(p$phi1,p$phi2),plot=TRUE,lag.max=40)
# How? Check theoretical acf again, now drawn as a curve with integer lags added as points:
z2 <- ar2.acf(a1=p$a1, a2=p$a2, plot=TRUE, subplot=TRUE)
points(as.numeric(names(z1$acf)), z1$acf, col=2)
abline(v=k, lty=2, col=2) # add vertical line at true peak
# Since the acf is evaluated at integer lags, it can appear that the dominant peak is located at 2*k for high ipv.


# BIFURCATIONS of the dominant peak of the acf
# Set up a grid with parameters covering a1/a2-space under the parabola
phi1 <- seq(-1.995,1.995,0.01)
phi2 <- seq(-0.995,0,0.01)
p <- expand.grid(phi1=phi1,phi2=phi2)

f <- function(phi1,phi2) phi1^2 + 4*phi2
z <- f(phi1=p[,1], phi2=p[,2])
inds <- which(z>0)
p <- p[-inds,] # remove values above the parabola

# Plot parameter space
ar2.plot(k=NULL, v=NULL)
points(p[,1], p[,2], col=2, cex=0.5, pch=".")

# Extract dominant period of theoretical acf derived from AR(2) coefficients.
# The lag.max options is set to 2*20, unless estimated quasi-period > 20, when lag.max is set to 2*length of quasi-period.
res <- numeric(nrow(p))
# Use ar2.period2, since we have generated phi1 = (1+a1) & phi2=a2
ar2.per <- ar2.period2(phi1=p[,1],phi2=p[,2],method="polyroot",output=list)$period

system.time({test <- for (i in 1:nrow(p)) {
	phi <- as.numeric(p[i,])
	lag.max <- ifelse(ar2.per[i] > 20, 2*ar2.per[i], 2*20)
	res[i] <- ar.acf(ar=as.numeric(phi), plot=FALSE, lag.max=lag.max)$'Dominant period'
}
})

dat <- cbind(p, acf.period=res, ar2.period=ar2.per)

# Create matrix for plotting.
mat <- tapply(dat[,3],dat[,1:2],c) # Notice clever usage of tapply
# Creata a list for image and contour plots:
g <- list(x = as.numeric(rownames(mat)), y = as.numeric(colnames(mat)), z = mat)

# Various plots
ar2.plot(k=NULL,v=NULL)
contour(g,levels=2:20,add=TRUE)

ar2.plot(k=seq(2,7,0.5),v=NULL, cex.axis=0.8, cex.contour=0.8, text.lab=FALSE)
contour(g,levels=2:20,add=TRUE)

ar2.plot(k=seq(2,7,0.5),v=NULL, cex.axis=0.8, cex.contour=0.8, text.lab=FALSE, ylim=c(-1,0), triangle=FALSE)
contour(g,levels=2:20,add=TRUE)

ar2.plot(k=NULL,v=NULL)
contour(g,levels=c(7,14),add=TRUE)

ar2.plot(k=c(3.5, 7, 14),ylim=c(-1,0),triangle=FALSE,width=8,height=5,text.lab=FALSE)
with(z <- subset(dat, dat[,3]==7|dat[,3]==14), points(z[,1], z[,2], col=2, pch="."))
with(z <- subset(dat, dat[,3]==3), points(z[,1], z[,2], col=3, pch="."))
with(z <- subset(dat, dat[,3]==4), points(z[,1], z[,2], col=4, pch="."))

ar2.plot(xlim=c(-.8,-.2),ylim=c(-1,-.75),text.lab=FALSE,cex.axis=0.8,cex.contours=0.8)
text(dat[,1],dat[,2],dat[,3],cex=0.5, col=ifelse(dat[,3]>7,2,ifelse(dat[,3]==7,4,1)))

ar2.plot(xlim=c(-2,-1.2),ylim=c(-1,-.85),text.lab=FALSE,cex.axis=0.8,cex.contours=0.8)
text(dat[,1],dat[,2],dat[,3],cex=0.5, col=ifelse(dat[,3]>7,2,ifelse(dat[,3]==7,4,1)))

filled.contour(g, levels=2:11, color.palette=terrain.colors, bty="l", main="Dominant ACF lags")

# UPDATE: How does all this compare with spectral and wavelet analysis???
# A: Spectral density does not show spurious double-period peaks in the critical regions, neither do wavelets.
