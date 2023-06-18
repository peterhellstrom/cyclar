library(cyclar)

# Find AR(2)-parameters
k <- c(3, 4, 5, 6, 7, 8) # Periodicity
v <- c(2, 2, 2, 3, 3, 1.5) # Variance

# Two different (but equivalent) functions for parameters:
ar2.parms(k = k, v = v, plot = TRUE, list = FALSE)
ar2.parms.old(k = k, v = v, plot = TRUE, output = cbind)

# ar2.parms.old has an extra argument, kmethod = c("royama","jenkins").
# Only "royama" is used in ar2.parms.
ar2.parms.old(k = k, v = v, kmethod = "jenkins", plot = TRUE)

# Add points to existing plot
ar2.plot(k = NULL, v = NULL)
pr <- ar2.parms(k = k, v = v)
points(1 + pr$a1, pr$a2, col = 2, cex = 2)

# Random parameter values
ar2.k.rand(n = 500, k = 6, spec = "norm")

# Test contour-functions
# Equally spaced values
ar2.k.gen(k = 6, length.out = 20, method = "x") # equal interval between x-values
ar2.k.gen(k = 6, length.out = 20, method = "x", from = 0.2, to = 0.5)

ar2.k.gen(k = 6, length.out = 20, method = "y") # equal interval between y-values
ar2.k.gen(k = 3, length.out = 30, method = "y", from = -1, to = -.5)

p <- ar2.k.gen(k = 3, length.out = 20, method = "x", from = -1, to = 0)
p <- ar2.k.gen(k = 3, length.out = 20, method = "x", from = -0.8, to = -0.2)
p <- ar2.k.gen(k = 3, length.out = 20, method = "x", from = -1)
p <- ar2.k.gen(k = 3, length.out = 20, method = "x", to = -0.8)
p <- ar2.k.gen(k = 3, length.out = 20, method = "x", from = 1) # error

p <- ar2.k.gen(k = 2, length.out = 20, method = "x", from = -1, to = 0)
p <- ar2.k.gen(k = 2, length.out = 20, method = "x", from = -0.8, to = -0.2)
p <- ar2.k.gen(k = 2, length.out = 20, method = "x", from = -1)
p <- ar2.k.gen(k = 2, length.out = 20, method = "x", to = -0.8)

p <- ar2.k.gen(k = 6, length.out = 20, method = "x", from = 0, to = 1)
p <- ar2.k.gen(k = 6, length.out = 20, method = "x", from = 0.1, to = 0.4)
p <- ar2.k.gen(k = 6, length.out = 20, method = "x", from = 0.55)
p <- ar2.k.gen(k = 6, length.out = 20, method = "x", to = 0.75)

# For variance (only possible to generate value on y-axis)
# Random values
ar2.v.rand(n = 500, v = 6, spec = "unif", min = -0.95, max = -0.55)

# Equal interval
ar2.v.gen(v = 6, length.out = 100, from = -.8, to = -.5)
ar2.v.gen(v = 6, length.out = 100, from = -.8, to = 0) # automatically adjusted, since to = 0 is too large
ar2.v.gen(v = 6, length.out = 100, to = -.5)
ar2.v.gen(v = 6, length.out = 100, from = -.25)

p <- ar2.parms(k = 8, v = 1.5)
ar2.period(p$a1, p$a2)
ar2.period2(p$phi1, p$phi2)


#### Test ar2.plot-function ####
ar2.plot()
ar2.plot(k = NULL, v = NULL)
ar2.plot(kcol.lab = 2)
ar2.plot(kcol.lab = 2, vcol.lab = 4, vlab.both = TRUE)
ar2.plot(triangle = FALSE, ylim = c(-1,0))
ar2.plot(text.lab = FALSE, par.name = "Omega")
ar2.plot(text.lab = FALSE, par.name = "Omega", xlab = expression(Omega[1]))
ar2.plot(text.lab = FALSE, triangle = FALSE, ylim = c(-1, 0), vlab.both = TRUE, kcol = "darkgrey", vcol = "lightgrey")
ar2.plot(text.lab = FALSE, formula = TRUE, formula.pos = c(-1, 1),par.name = "alpha")
ar2.plot(k = 5)
ar2.plot(k = 5, vlab = FALSE)
ar2.plot(k = 3:5, vlab = FALSE)
ar2.plot(k = 5, vlab = FALSE, triangle = FALSE)
ar2.plot(k = 5, vlab = FALSE, vcontours = FALSE)
ar2.plot(k = 5, v = 7, vlab = TRUE, vlab.both = TRUE, vregion = "triangle")
ar2.plot(k = 5, v = 7, kcol = 8, vlty = 1, vcol = "steelblue", vlab = TRUE, vlab.both = TRUE, vregion = "triangle")
ar2.plot(k = 5, v = 7, kcol = 8, vlty = 1, vcol = "steelblue", vlab = TRUE, vlab.both = TRUE, vregion = "triangle", cex.axis = 1.2)
ar2.plot(k = c(3,5), v = 7, kcol = 8, vlty = 1, vcol = "steelblue", vlab = TRUE, vlab.both = TRUE, vregion = "triangle", cex.axis = 1.2, cex.contours = 2)

# The Pink Floyd "Dark Side of the Moon" version of the AR(2)-plot
p <- ar2.parms(k = 8, v = 5, plot = FALSE)
ar2.period(p$a1,p$a2)
#png(filename = "DarkSide.png", res = 800, units = "in", width = 8, height = 4, pointsize = 11)
ar2.dsotm(a1 = p$a1, a2 = p$a2, n = 200, main = "Ticking away the moments that make up a dull day")

#### TEST contour functions ####
# Five functions for each kontour (k & v)
# ar2.*.add
# ar2.*.contour
# ar2.*.gen
# ar2.*.inv
# ar2.*

#### v-contours ####

# Example 1
ar2.v.gen(v = 4, plot = TRUE, by = 0.1)

# Example 2
ar2.plot(v = 6, vlab = TRUE, vlab.pos = 3, vcontours = FALSE)
ar2.v.add(v = 6, col = 2, lty = 1, lwd = 2)
ar2.k.add(k = 6, col = 4, lty = 1, lwd = 2)

# Example 3
z <- ar2.v.rand(v = 5, n = 100, spec = "norm", mean = -0.7, sd = 0.2)
points(1+z$a1, z$a2, pch = "+", col = 2)

# Example 4
v <- c(1.05, 1.25, 2, 3, 6, 10, 20)
ar2.plot(v = v, cex.contours = 0.75, vcol = "steelblue")

# Example 5
p1 <- ar2.v.gen(v = 2, vregion = "parabola", length.out = 50, plot = FALSE)
p2 <- ar2.v.gen(v = 5, vregion = "parabola", length.out = 50, plot = FALSE)
ar2.plot(k = NULL,v = c(2,5))
points(p1$phi1,p1$phi2)
points(p2$phi1,p2$phi2,col = 2)

# Example 6
ar2.v.rand(v = 2, vregion = "triangle", spec = "unif", min = -1, max = 1, n = 100)
ar2.v.rand(v = 5, vregion = "triangle", spec = "unif", min = -1, max = 1, n = 100)

# Example 7
# ar2.parms & ar2.period can handle vectors

k <- c(2,3,5,6,7)
v <- c(2,5,3,4,6)
ar2.plot(k = k, v = v, vregion = "triangle", vcol = "steelblue", vlty = 2, vlab.both = TRUE, formula = TRUE)
p <- ar2.parms(k = k, v = v)
p
points(1+p$a1, p$a2, pch = 16,col = 2)

ar2.period(a1 = p$a1, a2 = p$a2)

# Example 8
ar2.plot(v = 2,kcontours = FALSE,klab = FALSE,vcol = "steelblue")

# Example 9
ar2.plot(k = NULL, v = 6, vlab = FALSE, vlab.pos = 3, vcontours = FALSE)
z <- ar2.v.add(v = 4, vregion = "triangle", col = 3, lty = 3, lwd = 2)
points(1+z$a1, z$a2, pch = ".", col = 2)

#### k - contours ####

# Example 1
ar2.k.rand(k = 3, spec = "norm")
ar2.k.rand(k = 3, spec = "norm", method = "x")

# Example 2
# Example 1
ar2.k.rand(k = 3, spec = "norm", mean = -0.5, sd = 0.2)

# Example 3 (add variance contours)
p <- ar2.k.gen(k = 10, length.out = 20, from = -0.95)
p2 <- ar2.period(p$a1, p$a2, output = list)
ar2.plot(k = NULL, v = NULL)
for (i in 1:length(p2$ipv)) ar2.v.add(v = p2$ipv[i], vregion = "triangle", lty = 2)

###############
p <- ar2.parms(k = 3:5, v = rep(2,3), output = list)
ar2.period(p$a1, p$a2)

# Values from table 1 in Henden et al (2008) j appl ecol, p. 1089
ar2.period(a1 = c(-0.767955, 0, 0.5060442)-1, a2 = c(-0.589755, -0.707112, -0.67043))

ar2.period2(phi1 = c(-0.767955, 0, 0.5060442), phi2 = c(-0.589755, -0.707112, -0.67043))

#### Symmetry at -a1 and a1 ####
# Calculate period lengths for -a1 & a1

f <- function(a1,k) (-0.25*(1+a1)^2) * (1 + tan(2*pi/k)^2)
k <- function (a1,a2) (2*pi) / (acos((1+a1)/(2*sqrt(-a2))))

k <- 7
x <- seq(-2,2,0.1)-1
y <- f(a1 = x,k = k)

ar2.plot(k = k,v = NULL)
points(1+x,y)

z <- ar2.period(a1 = x, a2 = y)

ar2.plot.simple(a1 = z[,1],a2 = z[,2],ylim = c(f(a1 = x[1],k = 7),1))
curve(f(x-1,k = 7),add = T,col = 2,lty = 2)

###
k <- c(2,6,8,9) # quasi-period
v <- 6 # intrinsic process variance

ar2.plot(k = NULL,v = 6,kcontours = FALSE,klab = FALSE)

for (i in 1:length(k)) {

	p <- ar2.parms(k = k[i], v = v)
	a1 <- p$a1
	a2 <- p$a2
	curve(f(a1 = x-1,k = k[i]),col = i,add = T)
	points(1+a1,a2,pch = 16,col = 2)
	points(-(1+a1),a2,pch = 16,col = 2)

	a1.pos <- ar2.period(a1,a2)
	a1.neg <- ar2.period(-a1-2,a2)

	text(x = 2*cos(2*pi/a1.pos[,"period"]), y = -1, labels = round(a1.pos[,"period"],2), font = 2, cex = 1,col = 2)
	text(x = 2*cos(2*pi/a1.neg[,"period"]), y = -1, labels = round(a1.neg[,"period"],2), font = 2, cex = 1,col = 4)
}


#### Generate "uniformly" distributed a1 and a2 values ####
n <- 10000
out <- matrix(ncol = 2, nrow = n)
colnames(out) <- c("a1", "a2")

for (i in 1:n) {
	x <- runif(1, -2, 2) # (1+a1)
	a2.max <- -0.25*x^2
	a2 <- runif(1, -1, a2.max)
	out[i,] <- c(a1 = x-1,a2 = a2)
}

ar2.plot.simple()
points(1+out[,"a1"],out[,"a2"],pch = 16,col = 2)

par(mfrow = c(1,2))
	hist(1+out[,"a1"], col = "steelblue", breaks = 30, main = "(1+a1)")
	hist(out[,"a2"], col = "steelblue", breaks = 30, main = "a2")
par(mfrow = c(1,1))
