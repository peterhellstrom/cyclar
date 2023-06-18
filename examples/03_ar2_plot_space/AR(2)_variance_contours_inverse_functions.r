library(cyclar)

# The functions in this document are available in cyclar.

# The functions ar2.v.inv & ar2.k.inv were developed for plotting period - and variance contours in (1+a1),a2 space.
# The function for period(k)-contours can be readily found in Royama, and with some algebra it is easy to develop
# a formula for a2-values as a function of (1+a1) for k-contours. However, I haven't found any mention
# of the formula for variance contours. For variance contours, it appears that the inverse function [(1+a1) on the y-axis
# and a2 on the x-axis] has a fairly easy solution. I used this fact for e.g. ar2.plot. The variance contours are generated
# from the inverse function, and then the axes are shifted. The procedure is shown below.

# Note that ar2.v.inv is used in the function ar2.v.add.
# ar2.k.inv is included for completeness, but not yet used for any purpose in cyclar.

# Inverse functions for the k- and v-contours.
# x = a2, y = (1+a1)

# Notes on inverse functions:
# 1) replace f(x) with y
# 2) switch x's and y's
# 3) solve for y
# 4) replace y with f^-1(x)

ar2.v.inv <- function(x, v, na.rm=TRUE, out="y") {

	if (any(out %in% c("y","both")) == FALSE) stop("out must be one of = c('both','y')")

	#y2 <- (-1 + x) * (1 + v*(-1 + x^2)) / (v*(1+x))
	y2 <- x^2 - 2*x + 1 - ((1 - x )/((1 + x)*v))

	inds <- which(sign(y2)==-1)

	if (na.rm == FALSE) y <- sqrt(y2)
	else if (na.rm == TRUE) y <- sqrt(abs(y2))

	if (out == "both") list(y2=y2, y=y, NAs=inds, y2NAs=y2[inds])
	else if (out == "y") y
}

ar2.k.inv <- function(x, k) {
	# if k < 4, cos(2*pi/k) is negative
	# if k > 4, cos(2*pi/k) is positive
	# if k = 4, cos(2*pi/k is zero)
	#y <- -2*cos(2*pi/k)*sqrt(-x)
	y <- 2*cos(2*pi/k)*sqrt(-x) # The second solution
	y
}


# Plot the inverse functions
# First, set the (intrinsic process) variance
v <- 2

# Then creat plot with all the necessary components
plot(x=c(-1,1), y=c(-2,2), type="n", xlab=expression(beta[2]), ylab=expression(1 + beta[1]), main="Inverse functions")

lines(x=c(-1,1,-1,-1), y=c(2,0,-2,2),col=2) # Triangle
abline(h=0, lty=2)
curve(2*sqrt(-x), from=-1, to=0, n=1001, add=TRUE, col=2) # Parabola, positive values
curve(-2*sqrt(-x), from=-1, to=0, n=1001, add=TRUE, col=2) # Parabola, negative values
# Add v-contours for v
v.pos <- curve(ar2.v.inv(x=x, v=v, na.rm=TRUE), from=-sqrt(1-(1/v)), to=sqrt(1-(1/v)), col=4, lty=2, n=1001, add=TRUE)
v.neg <- curve(-ar2.v.inv(x=x, v=v, na.rm=TRUE), from=-sqrt(1-(1/v)), to=sqrt(1-(1/v)), col=4, lty=2, n=1001, add=TRUE)
# Add two k-contours, k=c(3,5)
k3 <- curve(ar2.k.inv(x=x, k=3), col=2, lty=2, n=1001, from=-1, to=0, add=TRUE)
k5 <- curve(ar2.k.inv(x=x, k=5), col=4, lty=2, n=1001, from=-1, to=0, add=TRUE)
# Add points at intersection between parabola and v-contours
# x.max is the a2-value at the intersection
# Plug x.max into ar2.v.inv to get the (1+a1)-value.
x.max <- (sqrt((27*v+1)/v)/(3^(3/2)*v)+1/v)^(1/3)-1/(3*v*(sqrt((27*v+1)/v)/(3^(3/2)*v)+1/v)^(1/3))-1
points(rep(x.max,2), c(ar2.v.inv(x=x.max,v=v),-ar2.v.inv(x=x.max,v=v)), pch=16, col=2)
# Add legend
legend("topright", legend=bquote(v == .(v)), col=4, lty=2, bty="n")

# ar2.v.inv has an na.rm option. Why?
# For x-values that produce y-values close to zero, the y2-values in ar2.v.inv is sometimes essentially zero
# but negative, i.e. -1.05e-9. Of course the sqrt-function returns an error (NaNs produced) in such cases. The na.rm=TRUE option
# squares the absolute values of y2. Check the two examples below (the out=c("y","both") option controls the output, set to both
# if you want both the y2 and y values to be printed, necessary in this case).

v <- 7
x1 <- seq(from=-sqrt(1-(1/v)), to=sqrt(1-(1/v)), length.out=1001)
y1 <- ar2.v.inv(x=x1,v=v,na.rm=FALSE,out="both")
str(y1)
sqrt(1-(1/v)); x1[y1$NAs]; y1$y2[y1$NAs]; y1$y[y1$NAs]

y2 <- ar2.v.inv(x=x1,v=v,na.rm=TRUE,out="both")
str(y2)
sqrt(1-(1/v)); x1[y2$NAs]; y2$y2[y2$NAs]; y2$y[y2$NAs]


# The next step is to construct a function that plots a v-contour in
# the "correct" space, that is x=(1+a1) & y=a2. This function is named ar2.v.add
# and is called e.g. by ar2.plot.

# Note on max- and min values of the function:
# The function is not defined outside +/- sqrt(1-(1/v))

# For vregion="triangle"
# Find the limits by taking the expression for y in ar2.v.inv, set equal to zero, and solve for x. Simplify.

# Mathematica code:
# vinv := sqrt((-1 + x) * (1 + v*(-1 + x^2)) / (v*(1+x)));
# z := Solve[vinv == 0,x];
# Simplify[z]

# Find y.max in parabola-region.
# For vregion="parabola", find the intersection between the v-contour and the parabola.
# The v-contours and the parabola intersect at this point
# Set the inverse functions for the contours and parabola equal to each other:
# sqrt((-1+x)*(1+v*(-1+x^2))/(v*(1+x))) = 2*sqrt(-x)
# Square both sides, and solve for x:
# (-1+x)*(1+v*(-1+x^2))/(v*(1+x)) = -4*x
# Three solutions, this is the one without imaginary parts:
# x = (sqrt((27*v+1)/v)/(3^(3/2)*v)+1/v)^(1/3)-1/(3*v*(sqrt((27*v+1)/v)/(3^(3/2)*v)+1/v)^(1/3))-1

# Mathamatica code:
#z := Solve[(-1+x)*(1+v*(-1+x^2))/(v*(1+x)) == -4*x, x];
#z
#z[[1]]

ar2.v.add <- function(v, vregion="parabola", n=1001, ...) {
	# y = a2, x = (1+a1)
	if (vregion == "triangle") {
		y.min <- -sqrt(1-(1/v))
		y.max <- -y.min
		yv <- seq(from=y.min, to=y.max, length.out=n)
	}
	if (vregion=="parabola") {
		y.min <- -sqrt(1-(1/v))
		y.max <- (sqrt((27*v+1)/v)/(3^(3/2)*v)+1/v)^(1/3)-1/(3*v*(sqrt((27*v+1)/v)/(3^(3/2)*v)+1/v)^(1/3))-1
		yv <- seq(from=y.min, to=y.max, length.out=n)
	}

	xv <- ar2.v.inv(x=yv,v=v,na.rm=TRUE)
	# xv only contains positive (1+a1)-values, but the curve is symmetric.
	points(xv,yv,type="l",...)
	points(-xv,yv,type="l",...) # plot negative (1+a1) values.
	invisible(list(a1=c(-xv-1,xv-1), a2=c(yv,yv), x=c(-xv,xv), y=c(yv,yv)))
}

# Test the function ar2.v.add (& ar.k.add)

# Example 1
ar2.plot.simple()
yv <- ar2.v.add(v=v, col=4, lty=1 ,lwd=1)
str(yv)
yk <- ar2.k.add(k=5, col=2, lty=2)
str(yk)

# Example 2
ar2.plot.simple()
#lines(y$x, y$y)
points(yv$x, yv$y, col=2, pch=".")

# Example 3, plot only in parabola-region
ar2.plot.simple()
ar2.v.add(v=1.1, col=2, lty=2)
ar2.v.add(v=2, col=2, lty=2)
ar2.v.add(v=2.1, col=4, lty=3)
ar2.v.add(v=4, col=4, lty=3)
ar2.v.add(v=7, col=4, lty=3)

# Example 4, plot in "full" triangle region
ar2.plot.simple()
ar2.v.add(v=1.1, vregion="triangle", col=2, lty=2)
ar2.v.add(v=2, vregion="triangle", col=2, lty=2)
ar2.v.add(v=2.1, vregion="triangle", col=4, lty=3)
ar2.v.add(v=4, vregion="triangle", col=4, lty=3)
ar2.v.add(v=7, vregion="triangle", col=4, lty=3)
ar2.k.add(k=5,col=2,lty=2)

# How did I get the expression for y in the ar2.v.inv function?

# The expression for the intrinsic process variance [Royama-style notation, e.g. phi1=(1+a1) and phi2=a2]
# (1-a2) / ((1-(1+a1)-a2)*(1-a2+(1+a1))*(1+a2)) = ipv
# Set a2 = y
# Set (1+a1) = x
# Set ipv = v

# Solve the following equation for x
# (1-y) / ((1-x-y)*(1-y+x)*(1+y)) = v

# Multiply & simplify the following part in the denominator: (1-x-y)*(1-y-+x)=
# y^2 - 2*y - x^2 + 1

# The equation is now:
# (1 - y) / ((y^2 - 2*y - x^2 + 1) * (1 + y)) = v

# Multiply both sides with (y^2 - 2*y - x^2 + 1)
# (1 - y)/(1 + y) = v*y^2 - 2*v*y - v*x^2 + v

# Divide both sides with v:
# (1 - y) / ((1 + y)*v) = y^2 - 2*y - x^2 + 1

# Rearrange & isolate x^2
# x^2 = y^2 - 2*y + 1 - ((1 - y )/((1 + y)*v))

# The solution is:
# x = +/- sqrt(y^2 - 2*y + 1 - ((1 - y )/((1 + y)*v)))

# Note that this gives us the (1+a1) values, if v & a2 is known.

# Solve the above equation for y, and you have the equation
# for the v-contours as a function of x (1+a1). Solutions exists,
# but are difficult...

# Switch x's and y's, and we have the formula that is found in ar2.v.inv.
