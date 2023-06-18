################################################################################
library(cyclar)
library(ellipse)
################################################################################
# Simulate correlated variables
# This approach is probably very similar to the well-known Moran effect

k <- 5
v <- 3
p <- ar2.parms(k=k,v=v)
mu <- 2.1
a1 <- p$a1
a2 <- p$a2
ar2.var <- 0.2

n <- 1000
#npatch <- 1000 # For simulations of method accuracy
npatch <- 500
n.start <- 1000

# Try changing rho
# The both methods should give very similar results if rho=0,
# in which case there are no correlation between patches (diagonal covariance matrix)
rho <- 0.8
#rho <- 0

xpatch <- ar2.sim.npatch(n=n, npatch=npatch, rho=rho, mu=mu, a1=a1, a2=a2, sd=sqrt(ar2.var), n.start=n.start, plot=FALSE)
xsim <- sapply(1:npatch, function(i) ar2.sim(n=n, mu=mu, a1=a1, a2=a2, sd=sqrt(ar2.var), n.start=n.start))

xpatch.stat <- t(sapply(1:npatch, function(i) ts.summary(xpatch$x[,i], plot=FALSE)))
xsim.stat <- t(sapply(1:npatch, function(i) ts.summary(xsim[,i], plot=FALSE)))

colMeans(xpatch.stat)
colMeans(xsim.stat)

# Expected sd
sqrt(ar2.ipv(a1,a2,sd=sqrt(ar2.var)))

# Check correlation between time series (do this only if you have simulated few patches, i.e. npatch is small)
#pairs(xpatch, panel=panel.smooth)
#cor(xpatch$ts.x)
#cov(xpatch$ts.x)
#var(xpatch$ts.x)

# Estimate AR(2)-coefficients from simulated series
xpatch.coef <- t(sapply(1:npatch, function(i) ar.yw(x=xpatch$x[,i], order.max=2)$ar))
xsim.coef <- t(sapply(1:npatch, function(i) ar.yw(xsim[,i], order.max=2)$ar))

# Compare expected and observed values of simulations
ar2.plot(ylim=c(-1,0),k=k,v=v, triangle=FALSE)
points(xsim.coef[,1], xsim.coef[,2], pch="+", col=3)
points(xpatch.coef[,1], xpatch.coef[,2], pch="+", col=4)
points(1+a1,a2,col=2,pch=16,cex=1.3) # Simulation input
legend("topright",c("uncorrelated","correlated"), pch="+", col=3:4, bty="n")

# Check distribution of estimated AR(2)-coefficients
par(mfrow=c(2,1))
plot(density(xpatch.coef[,1]),main="Direct density-dependence")
lines(density(xsim.coef[,1]),col=2)
legend("topleft",c("uncorrelated","correlated"), col=c(2,1), lty=1, bty="n")

plot(density(xpatch.coef[,2]),main="Delayed density-dependence")
lines(density(xsim.coef[,2]),col=2)
legend("topleft",c("uncorrelated","correlated"), col=c(2,1), lty=1, bty="n")
par(mfrow=c(2,1))

# dput(ts.intersect(xpatch[,1:2]),"wavelet_coherency_example")


# Simulate two different habitats ----
# Remember that rho can take any value between -1 & 1!

z2 <- ar2.sim.npatch(n=100, npatch=2, rho=0.75, mu=0, a1=0, a2=-.7, sd=c(0.2,0.5), n.start=500, plot=TRUE)

# Perfect positive correlation
z2a <- ar2.sim.npatch(n=100, npatch=2, rho=1, mu=0, a1=0, a2=-.7, sd=c(0.2,0.5), n.start=500, plot=TRUE)
cor(z2a$x)
cor(z2a$w)
# Perfect negative correlation
z2b <- ar2.sim.npatch(n=100, npatch=2, rho=-1, mu=0, a1=0, a2=-.7, sd=c(0.2,0.5), n.start=500, plot=TRUE)
cor(z2b$x)
cor(z2b$w)
# Independent
z2c <- ar2.sim.npatch(n=100, npatch=2, rho=0, mu=0, a1=0, a2=-.7, sd=c(0.2,0.2), n.start=500, plot=TRUE)
cor(z2c$x)
cor(z2c$w)

# Let mu vary between habitats
z3 <- ar2.sim.npatch(n=500, npatch=2, rho=0.8, mu=c(0.5,2), a1=0, a2=-.7, sd=c(0.2,0.5), n.start=500, plot=TRUE)
cor(z3$x)
cov(z3$x)

# Plot correlation between habitats
plot(as.numeric(z3$x[,1]),as.numeric(z3$x[,2]), xlab="Habitat 1", ylab="Habitat 2", font.lab=2)
lines(ellipse(cov(z3$x), centre=c(0.5,2), level=0.90), lty=2) # add ellipse for confidence region, from package ellipse

# Cross-correlation between habitats
ccf(z3$x[,1],z3$x[,2], main="Cross-correlation between habitat 1 & 2")

# Add "disturbance"
# Input mu as vector
# Time-varying?
