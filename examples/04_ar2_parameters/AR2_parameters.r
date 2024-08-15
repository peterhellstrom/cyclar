library(cyclar)

# 1) Find parameters (1 + a1) and a2, given cycle period (k) and intrinsic process variance (v) ----

k <- c(2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10) # Cycle period
v <- c(1.05, 1.1, 1.25, 1.5, 2, 3, 4, 6, 8, 10) # Variance (intrinsic process variance)

## Use the function ar2.parms ----
inp <- expand.grid(k=k, v=v)
nrow(inp)

## Find parameters in AR2-space ----
p <- ar2.parms(k=inp$k, v=inp$v, output=list)
# Check "backward compatibility"
ar2.period(a1=p$a1, a2=p$a2)

## Add parameters to plot ----
ar2.plot(k=k, v=v, klab=FALSE, vlab=FALSE) # No contour labels
# Add points for intersections between period and variance contours:
points(1+p$a1,p$a2, pch=16, col=2, cex=1.2)

# 2) Plot triangle-space and add some points ----
a1 <- c(-0.767955, 0, 0.5060442) - 1
a2 <- c(-0.589755, -0.707112, -0.67043)

k <- c(2.5, 3, 3.5, 4, 5, 6, 7, 8, 9) # Period lengths
v <- c(2,4,10) # intrinsic process variance

ar2.plot(k=k,v=v)
points(1+a1,a2,col=2,pch=16,cex=1.3)

ar2.plot()
points(1+a1,a2,col=2,pch=16,cex=1.3)

ar2.plot(k=k,vcontours=FALSE,vlab=FALSE)
points(1+a1,a2,col=2,pch=16,cex=1.3)

ar2.plot(kcontours=FALSE, vcontours=FALSE)
points(1+a1,a2,col=2,pch=16,cex=1.3)

ar2.plot(kcontours=FALSE, vcontours=TRUE)
points(1+a1,a2,col=2,pch=16,cex=1.3)

# 3) Plot only parts of the triangle space (specify xlims and ylims) ----
ar2.plot(klab=FALSE, vlab=FALSE, xlim=c(0.075,2),ylim=c(-1,0), triangle=FALSE)

# 4) Find contour-lines and points ----

ar2.k.gen(k=4, method="y", length.out=50)
ar2.k.gen(k=8, method="y", length.out=50)

c35 <- ar2.k.gen(k=3.5, method="y", length.out=30)
ar2.period(a1=c35$a1, a2=c35$a2)

# 5) Mean, median, mu and a0 ----
# sigma (sd) vs. intrinsic process variance

# Generate parameters for 3,4,5-year cycles, all with intrinsic process variance=2 ----
inp <- list(k = 3:5, v = rep(2,3))
p <- ar2.parms(k=inp$k, v=inp$v, output=list)
p

# Check calculations (calculate period & ipv from AR(2) parameters ----
ar2.period(a1 = p$a1, a2=p$a2, method="polyroot")
ar2.period(a1 = p$a1, a2=p$a2, method="eigen")
ar2.period2(phi1 = p$phi1, phi2=p$phi2, method="eigen")

# Assume that mean on log-scale is 2.1
# Find the corresponding intercept a0 in the recursive equation,
# i.e. x[t] = a0 + (1+a1)[t-1] + a2[t-2] OR x[t] = phi0 + phi1[t-1] + phi2[t-2]
a0 <- ar2.intercept(mu = rep(2.1, 3), a1 = p$a1, a2=p$a2)
phi0 <- ar2.intercept2(mu = rep(2.1, 3), phi1 = p$phi1, phi2=p$phi2)

## Check, convert back to check mu ----
ar2.intercept(a0 = a0, a1 = p$a1, a2=p$a2)
ar2.intercept2(phi0 = phi0, phi1 = p$phi1, phi2=p$phi2)

# Fred's skua paper
ar2.intercept(a0 = c(1.1, 1.6, 1.9), a1 = p$a1, a2=p$a2)
ar2.intercept(mu = c(1.1, 1.6, 1.9), a1 = p$a1, a2=p$a2) # which one is correct?

# TO DO:
# Simulate with ar2.sim, ar2.sim2(try different methods)

# Plot distribution of simulated values on both log and natural scales
# add mean, median to distributions

# Simulate over a range of variance values
# Plot variance on x-scale, and exp(x[t]) on y-scale
# Add expected and observed median and mean on both natural and log-scale

# Estimation of sigma/variance?
